use std::hash::{Hash, Hasher};

use super::*;
use rustc_hash::FxHasher;
use smallvec::SmallVec;

impl Mol {
    pub fn ecfp(&self, radius: u8) -> u64 {
        // however, benchmarking this shows no significant change using existing microbenches for
        // ecfp
        let mut implicit = self.clone();
        implicit.implicit_hs();

        // original atom identifiers (radius 0)
        let mut atom_ids = FxHashMap::<NodeIndex, u64>::default();
        let mut features = FxHashSet::<u64>::default();
        features.reserve(radius as usize * implicit.graph.node_count());

        // compute initial properties
        for atom_idx in implicit.graph.node_indices() {
            let mut hasher = FxHasher::default();
            // properties are the daylight atom invariants:
            // 1. the number of heavy neighbors
            implicit.graph.neighbors(atom_idx).count().hash(&mut hasher);
            // 2. valence - # of hydrogens
            implicit.total_bonds(atom_idx).hash(&mut hasher);
            // 3. atomic number
            implicit.graph[atom_idx].atomic_num().hash(&mut hasher);
            // 4. atomic mass
            let mass = implicit.graph[atom_idx].exact_mass();
            let int_mass: u64 = mass.to_bits();
            int_mass.hash(&mut hasher);
            // 5. atomic charge
            implicit.charge(&atom_idx).hash(&mut hasher);
            // 6. number of hydrogens
            implicit.h_count(&atom_idx).hash(&mut hasher);
            // 7. is in ring? -- NOTE: skipping for speed

            // store initial value
            atom_ids.insert(atom_idx, hasher.finish());
        }

        features.extend(atom_ids.values().cloned());

        let mut neighbors = SmallVec::<[(u64, u64); 4]>::default();
        for it in 1..=radius {
            let mut tmp_ids = FxHashMap::<NodeIndex, u64>::default();
            for atom_idx in implicit.graph.node_indices() {
                let mut hasher = FxHasher::default();
                for e in implicit
                    .graph
                    .edges(atom_idx)
                    .map(|edge| (*edge.weight() as u64, atom_ids[&edge.target()]))
                {
                    neighbors.push(e)
                }
                neighbors.sort_unstable();
                (it as u64).hash(&mut hasher);
                atom_ids[&atom_idx].hash(&mut hasher);
                for (weight, neigh_id) in neighbors.drain(..) {
                    weight.hash(&mut hasher);
                    neigh_id.hash(&mut hasher);
                }
                tmp_ids.insert(atom_idx, hasher.finish());
            }
            atom_ids = tmp_ids;
            features.extend(atom_ids.values().cloned());
        }

        // not removing duplicates because who cares

        let mut output = features.into_iter().collect::<Vec<_>>();
        output.sort_unstable();
        let mut hasher = FxHasher::default();
        for block in output {
            block.hash(&mut hasher);
        }

        hasher.finish()
    }

    /// Constructs a molecular hash compatible with InChIKeys in terms of storage from the
    /// [`Mol::ecfp`].
    pub fn ecfp_hash(&self, radius: u8) -> [u8; 14] {
        let ecfp = self.ecfp(radius).to_le_bytes();
        debug_assert_eq!(ecfp.len(), 8);
        let mut hash = [0u8; 14];

        // we want an easily representable hash
        // we are restricting ourselves to ASCII codes between ';' and 'Z'
        let mut stash = 0u32;
        for i in 0..8 {
            // if we zero out the first tree bits and add ';', then we're in our range
            // we will add these missing bits back in later
            hash[i] = (ecfp[i] & 0b0001_1111) + b';';

            // store those three bits for later use.
            // first create a u8 with 3 MSB in the LSB position
            // cast this as a u32, which 0-pads on the left
            // then shift over by 3*i to make sure we don't clobber bits that already exist
            stash |= (((ecfp[i] & 0b1110_0000) >> 5) as u32) << (3 * i);
        }

        // now add back in the three bits we were missing
        // we had 8*3 = 24 bits that are relevant
        // the first 20 bits we are going to store in the next 4 bytes
        for byte in 0..4 {
            // take 5 bits
            // we zero out LSB then mask out irrelevant parts
            let bits = stash >> (byte * 5) & 0b0001_1111;
            assert!(bits <= 31);
            hash[8 + byte as usize] = u8::try_from(bits).unwrap() + b';';
        }

        // now we just have 4 bits left to deal with, they'll get their own byte
        let bits = stash >> 20;
        assert!(bits <= 31);
        hash[12] = u8::try_from(bits).unwrap() + b';';

        // we have one byte left that we can store values from 0-31
        // carbon and hydrogen count are informative but probably won't fall in this range.
        // so we'll do nitrogen counts, no real reason but feels fine

        let n_count = std::cmp::min(
            self.graph.node_weights().filter(|&&x| x == Atom::N).count(),
            31,
        );
        hash[13] = n_count as u8 + b';';
        hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collision() {
        let m1 = MolBuilder::from_smiles("CNOC").unwrap().build();
        let m2 = MolBuilder::from_smiles("CONC").unwrap().build();
        assert_eq!(m1.ecfp_hash(2), m2.ecfp_hash(2));

        let smiles = std::fs::read_to_string("test_files/test_compounds.smi").unwrap();

        let mut hashes = FxHashMap::default();
        for smiles in smiles.lines() {
            let mol = MolBuilder::from_smiles(smiles).unwrap().build();
            let old = hashes.insert(mol.ecfp_hash(4), smiles);
            assert!(old.is_none(), "old: {}\ncurr: {}", old.unwrap(), smiles);
        }
    }

    #[test]
    fn explicit_h() {
        let pubchem_canonical = MolBuilder::from_smiles("CC(C)CC1C(=O)NC(C(=O)N2CCCC2C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N1)CCCN)C(C)C)CC3=CC=C(C=C3)O)CCC(=O)N)CC(=O)N)CC4=CC=CC=C4)CC5=CC=CC=C5)CC6=CC=CC=C6").unwrap().build();
        let pubchem_isomeric = MolBuilder::from_smiles("CC(C)C[C@H]1C(=O)N[C@@H](C(=O)N2CCC[C@H]2C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCCN)C(C)C)CC3=CC=C(C=C3)O)CCC(=O)N)CC(=O)N)CC4=CC=CC=C4)CC5=CC=CC=C5)CC6=CC=CC=C6").unwrap().build();
        let seq2nrp = MolBuilder::from_smiles("C=6C=CC=CC=6CC1NC(=O)C(CC(C)C)NC(=O)C(CCCN)NC(=O)C(C(C)C)NC(=O)C(CC5=CC=C(O)C=C5)NC(=O)C(CCC(N)=O)NC(=O)C(CC(N)=O)NC(=O)C(CC4=CC=CC=C4)NC(=O)C(CC3=CC=CC=C3)NC(=O)C2N(CCC2)C1=O").unwrap().build();

        assert_eq!(pubchem_canonical.ecfp_hash(8), seq2nrp.ecfp_hash(8));
        assert_eq!(
            pubchem_isomeric.ecfp_hash(8),
            pubchem_canonical.ecfp_hash(8)
        );
        assert_eq!(pubchem_isomeric.ecfp_hash(8), seq2nrp.ecfp_hash(8));
    }
}
