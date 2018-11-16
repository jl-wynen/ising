/**
 * Rust implementation of the Ising Model simulation.
 */

use std::io::prelude::*;
use std::fs;
use std::path::Path;
use std::ops::{Index, IndexMut};
use std::time::Instant;
use std::env;

extern crate rand;
use rand::prelude::*;


//--------------------------
// Set run parameters here.

const NTHERM_INIT: usize = 1000;  // number of thermalisation sweeps in the beginning
const NTHERM: usize = 1000;  // number of thermalisation sweeps per temperature
const NPROD: usize = 10000;  // number of production sweeps (with measurements) per temperature

const NX: usize = 4;  // number of lattice sites in x direction
const NY: usize = 3;  // number of lattice sites in y direction
const LATSIZE: usize = NX*NY;  // total lattice size

/// Return a vector of temperatures to run the simulation with.
fn list_temperatures() -> Vec<f64> {
    let mut temperatures: Vec<f64> = Vec::new();
    for i in 0..12 {
        temperatures.push((i as f64 + 1.)*0.5);
    }
    return temperatures
}

// End of run parameters.
//------------------------

/// Helper struct to handle a random number generator.
struct Rng {
    rng: StdRng,
}

impl Rng {
    /// Create an instance of Rng from a given seed.
    fn from_seed(seed: [u8; 32]) -> Rng {
        Rng{rng: StdRng::from_seed(seed)}
    }

    /// Generate a random index into a configuration.
    fn gen_index(&mut self) -> usize {
        use rand::Rng;
        self.rng.gen_range(0, LATSIZE)
    }

    /// Generate a random spin, one of {-1, +1}.
    fn gen_spin(&mut self) -> i32 {
        use rand::Rng;
        match self.rng.gen_range(0, 2) {
            0 => -1,
            _ => 1,  // 1 is the only other possibility
        }
    }

    /// Generate a random double in [0, 1].
    fn gen_real(&mut self) -> f64 {
        use rand::Rng;
        self.rng.gen_range(0., 1.)
    }
}

/// Hold a spin configuration on the lattice.
struct Configuration {
    /// The actual configuration, +1 for spin up, -1 for spin down.
    cfg: [i32; LATSIZE],

    /// List nearest neighbour indices for each site.
    /**
     * Neighbours for site i are stored at (4*i+0)...(4*i+3) in the order
     * x+1, x-1, y+1, y-1.
     */
    neighbours: [usize; 4*LATSIZE],
}

impl Configuration {
    /// Create a random configuration.
    fn random(rng: &mut Rng) -> Configuration {
        let mut cfg = Configuration{cfg: [0; LATSIZE],
                                    neighbours: make_neighbour_list()};

        for site in &mut cfg.cfg.iter_mut() {
            *site = rng.gen_spin();
        }

        cfg
    }
}

impl Index<usize> for Configuration {
    type Output = i32;

    /// Read spin at site idx.
    fn index(&self, idx: usize) -> &i32 {
        return &self.cfg[idx];
    }
}

impl IndexMut<usize> for Configuration {
    /// Modify spin at site idx.
    fn index_mut(&mut self, idx: usize) -> &mut i32 {
        return &mut self.cfg[idx];
    }
}

/// Store Monte-Carlo history of observables.
struct Observables {
    energy: Vec<f64>,
    magnetisation: Vec<f64>,
}

/// Return a list of nearest neighbour indices for use as neighbours in Configuration.
fn make_neighbour_list() -> [usize; 4*LATSIZE] {
    let mut indices: [usize; 4*LATSIZE] = [0; LATSIZE*4];

    for y in 0..NY {
        for x in 0..NX {
            indices[(y*NX+x)*4 + 0] = if x == NX-1 { y*NX } else { y*NX + x+1 };
            indices[(y*NX+x)*4 + 1] = if x == 0 { y*NX + NX-1 } else { y*NX + x-1 };
            indices[(y*NX+x)*4 + 2] = if y == NY-1 { x } else { (y+1)*NX + x };
            indices[(y*NX+x)*4 + 3] = if y == 0 { (NY-1)*NX + x } else { (y-1)*NX + x };
        }
    }

    indices
}

/// Create the output data directory and write the temperature file.
/**
 * Deletes the directory and all its contents if it exists.
 */
fn prepare_datadir(dirname: &Path, temperatures: &Vec<f64>) {
    if dirname.exists() {
        println!("Data directory '{}' exists, removing!", dirname.display());
        fs::remove_dir_all(dirname).unwrap();
    }
    fs::create_dir_all(dirname).unwrap();

    let mut tempfile = fs::File::create(dirname.join("temperatures.dat")).unwrap();
    for (i, temp) in temperatures.iter().enumerate() {
        write!(tempfile, "{}: {}\n", i, temp);
    }
}

/// Write observables to a data file.
fn write_observables(fname: &Path, obs: &Observables) {
    let mut obsfile = fs::File::create(fname).unwrap();

    for energy in obs.energy.iter() {
        write!(obsfile, "{} ", energy);
    }
    write!(obsfile, "\n");

    for magn in obs.magnetisation.iter() {
        write!(obsfile, "{} ", magn);
    }
    write!(obsfile, "\n");
}

/// Evaluate the Hamiltonian on a configuration.
fn hamiltonian(cfg: &Configuration) -> i32 {
    let mut energy: i32 = 0;

    for (idx, site) in cfg.cfg.iter().enumerate() {
        energy += site * (cfg[cfg.neighbours[4*idx]]
                          + cfg[cfg.neighbours[4*idx+1]]
                          + cfg[cfg.neighbours[4*idx+2]]
                          + cfg[cfg.neighbours[4*idx+3]]);
    }

    return -energy;
}

/// Compute the magnetisation on a configuration.
fn magnetisation(cfg: &Configuration) -> f64 {
    return cfg.cfg.iter().sum::<i32>() as f64 / LATSIZE as f64;
}

/// Compute the change in energy if the spin at site idx were flipped.
fn delta_e(cfg: &Configuration, idx: usize) -> i32 {
    return 2*cfg[idx] * (cfg[cfg.neighbours[4*idx]]
                         + cfg[cfg.neighbours[4*idx+1]]
                         + cfg[cfg.neighbours[4*idx+2]]
                         + cfg[cfg.neighbours[4*idx+3]]);
}

/// Evolve a configuration in Monte-Carlo time.
/**
 * Flips spins at random sites nsweep*NX*NY times and accepting or
 * rejecting the change using the Metropolis-Hastings algroithm.
 * Measures observables every NX*NY steps, i.e. once per sweep.
 *
 * cfg and energy must be set before calling the function.
 * Upon return, they contain the final configuration and energy.
 * Returns the number of accepted spin flips.
 */
fn evolve(cfg: &mut Configuration, energy: &mut f64, beta: f64,
          rng: &mut Rng, nsweep: usize, mut obs: Option<&mut Observables>) -> usize {
    // running number of accepted spin flips
    let mut naccept: usize = 0;

    for _sweep in 0..nsweep {
        for _step in 0..LATSIZE {
            let idx = rng.gen_index();  // flip spin at this site

            let delta = delta_e(&cfg, idx);  // proposed change in energy

            // Metropolis-Hastings accept-reject
            // The first check is not necessary for this to be correct but avoids
            // evaluating the costly exponential and RNG.
            if delta <= 0 || (-beta*(delta as f64)).exp() > rng.gen_real() {
                cfg[idx] *= -1;
                *energy += delta as f64;
                naccept += 1;
            }
            // else: discard
        }

        // measure observables if an instance of Observables is given.
        if let Some(o) = &mut obs {
            o.energy.push(*energy);
            o.magnetisation.push(magnetisation(&cfg));
        }
    }

    return naccept;
}

fn main() {
    // parse command line arguments
    let args: Vec<String> = env::args().collect();
    let datadir = if args.len() == 2 {
        Path::new(&args[1])
    }
    else {
        Path::new("./data")
    };

    // prepare output directory
    let temperatures = list_temperatures();
    prepare_datadir(&datadir, &temperatures);

    // one rng for all purposes
    let mut rng = Rng::from_seed([138; 32]);

    // initial condition (hot start)
    let mut cfg = Configuration::random(&mut rng);
    let mut energy = 0.0;  // does not matter for initial thermalisation

    // start measuring time, the above doesn't count
    let start_time = Instant::now();

    // initial thermalisation
    let naccept = evolve(&mut cfg, &mut energy, 1./temperatures[0], &mut rng, NTHERM_INIT, None);
    println!("Initial thermalisation acceptance rate: {}", (naccept as f64)/((NTHERM_INIT*LATSIZE) as f64));

    for (i, temp) in temperatures.iter().enumerate() {
        println!("Running for temperature {}", temp);
        let beta = 1./temp;
        energy = hamiltonian(&cfg) as f64;

        // re-thermalise
        let naccept = evolve(&mut cfg, &mut energy, beta, &mut rng, NTHERM, None);
        println!("  Thermalisation acceptance rate: {}", (naccept as f64)/((NTHERM*LATSIZE) as f64));

        // measure
        let mut obs = Observables{energy: Vec::new(), magnetisation: Vec::new()};
        let naccept = evolve(&mut cfg, &mut energy, beta, &mut rng, NPROD, Some(&mut obs));
        println!("  Production acceptance rate: {}", naccept as f64 / (NPROD*LATSIZE) as f64);

        write_observables(&datadir.join(format!("{}.dat", i)), &obs);
    }

    let duration = start_time.elapsed();
    println!("Duration in wall clock time: {}s", duration.as_secs() as f64
             + (0.001*duration.subsec_millis() as f64));
}
