use paths::*;

#[test]
fn parse_dirpath() {
    assert!("./fake_path".parse::<DirPath>().is_err());
    assert!("./src/".parse::<DirPath>().is_ok());
    assert!("../dereplicate".parse::<DirPath>().is_ok());
    assert!("Cargo.toml".parse::<DirPath>().is_err());
    assert!("../../target".parse::<DirPath>().is_ok());
    assert!("src".parse::<DirPath>().is_ok());
    assert!("./../".parse::<DirPath>().is_ok());
}

#[test]
fn parse_new() {
    assert!("./fake_path".parse::<NewDirPath>().is_ok());
    assert!("./src/".parse::<NewDirPath>().is_err());
    assert!("../dereplicate".parse::<NewDirPath>().is_err());
    assert!("Cargo.toml".parse::<NewDirPath>().is_err());
    assert!("../../target".parse::<NewDirPath>().is_err());
    assert!("src".parse::<NewDirPath>().is_err());
    assert!("./../".parse::<NewDirPath>().is_err());

    assert!("./src".parse::<NewFilePath>().is_err());
    assert!("./fake_path.rs".parse::<NewFilePath>().is_ok());
    assert!("../simple_app.rs".parse::<NewFilePath>().is_ok());
    assert!("./target".parse::<NewFilePath>().is_ok());
    assert!("./src/lib.rs".parse::<NewFilePath>().is_err());
    assert!("src/lib.rs".parse::<NewFilePath>().is_err());
    assert!("../../Cargo.toml".parse::<NewFilePath>().is_err());
    assert!("../constants/Cargo.toml".parse::<NewFilePath>().is_err());
}

#[test]
fn parse_filepath() {
    assert!("./src".parse::<FilePath>().is_err());
    assert!("./fake_path.rs".parse::<FilePath>().is_err());
    assert!("../simple_app.rs".parse::<FilePath>().is_err());
    assert!("./target".parse::<FilePath>().is_err());

    assert!("./src/lib.rs".parse::<FilePath>().is_ok());
    assert!("src/lib.rs".parse::<FilePath>().is_ok());
    assert!("../../Cargo.toml".parse::<FilePath>().is_ok());
    assert!("Cargo.toml".parse::<FilePath>().is_ok());
    assert!("../constants/Cargo.toml".parse::<FilePath>().is_ok());
}
