//! The run's log, one line per event as `[HH:MM:SS] message`, laid out as
//! ltr_library.py's: the screen (stderr) shows the stage headings, the results
//! and warnings; `-v` adds every step's detail line; `--log FILE` gets all of
//! them whatever the screen shows.

use std::io::Write;
use std::sync::{Mutex, OnceLock};

struct Log {
    verbose: bool,
    file: Option<Mutex<std::fs::File>>,
}

static LOG: OnceLock<Log> = OnceLock::new();

#[derive(PartialEq, Clone, Copy)]
pub enum Level {
    /// stage headings and results: always on the screen
    User,
    /// each step's numbers: on the screen with -v
    Detail,
    /// always on the screen
    Warn,
}

/// Set the screen level and the log file; before this (or without it) every line goes to the screen.
pub fn init(verbose: bool, file: Option<&str>) -> anyhow::Result<()> {
    let file = match file {
        Some(p) => Some(Mutex::new(
            std::fs::File::create(p).map_err(|e| anyhow::anyhow!("cannot create --log {}: {}", p, e))?,
        )),
        None => None,
    };
    LOG.set(Log { verbose, file }).map_err(|_| anyhow::anyhow!("the log was set up twice"))
}

/// The local wall-clock time as `[HH:MM:SS]`.
pub fn stamp() -> String {
    // SAFETY: time() and localtime_r() only write into the tm passed to them
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    unsafe {
        let t = libc::time(std::ptr::null_mut());
        libc::localtime_r(&t, &mut tm);
    }
    format!("[{:02}:{:02}:{:02}]", tm.tm_hour, tm.tm_min, tm.tm_sec)
}

pub fn write(level: Level, msg: &str) {
    let line = format!("{} {}", stamp(), msg);
    match LOG.get() {
        None => eprintln!("{}", line),
        Some(log) => {
            if level != Level::Detail || log.verbose {
                eprintln!("{}", line);
            }
            if let Some(f) = &log.file {
                if let Ok(mut f) = f.lock() {
                    let _ = writeln!(f, "{}", line);
                }
            }
        }
    }
}

/// A count with thousands separators: 18582025 -> 18,582,025.
pub fn n(x: usize) -> String {
    let s = x.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// A stage heading or a result line (always shown).
macro_rules! log_user { ($($t:tt)*) => { $crate::log::write($crate::log::Level::User, &format!($($t)*)) } }
/// A step's detail line (shown with -v; always in --log).
macro_rules! log_detail { ($($t:tt)*) => { $crate::log::write($crate::log::Level::Detail, &format!("  {}", format!($($t)*))) } }
/// A warning (always shown).
macro_rules! log_warn { ($($t:tt)*) => { $crate::log::write($crate::log::Level::Warn, &format!("  warning: {}", format!($($t)*))) } }
