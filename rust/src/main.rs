use std::io::prelude::*;
use std::fs;
use std::path::Path;
use std::ops::{Index, IndexMut};
use std::time::Instant;

extern crate rand;
use rand::prelude::*;

const NX: usize = 4;
const NY: usize = 3;
const LATSIZE: usize = NX*NY;

const NTHERM_INIT: usize = 1000;
const NTHERM: usize = 1000;
const NPROD: usize = 10000;

struct Rng {
    rng: StdRng,
}

impl Rng {
    fn from_seed(seed: [u8; 32]) -> Rng {
        Rng{rng: StdRng::from_seed(seed)}
    }

    fn gen_index(&mut self) -> usize {
        use rand::Rng;
        self.rng.gen_range(0, LATSIZE)
    }

    fn gen_spin(&mut self) -> i32 {
        use rand::Rng;
        match self.rng.gen_range(0, 2) {
            0 => -1,
            _ => 1,  // 1 is the only other possibility
        }
    }

    fn gen_real(&mut self) -> f64 {
        use rand::Rng;
        self.rng.gen_range(0., 1.)
    }
}

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

struct Configuration {
    cfg: [i32; LATSIZE],
    neighbours: [usize; 4*LATSIZE],
}

impl Configuration {
    fn random(rng: &mut Rng) -> Configuration {
        let mut cfg = Configuration{cfg: [0; LATSIZE],
                                    neighbours: make_neighbour_list()};

        for site in &mut cfg.cfg {
            *site = rng.gen_spin();
        }

        cfg
    }
}

impl Index<usize> for Configuration {
    type Output = i32;

    fn index(&self, idx: usize) -> &i32 {
        return &self.cfg[idx];
    }
}

impl IndexMut<usize> for Configuration {
    fn index_mut(&mut self, idx: usize) -> &mut i32 {
        return &mut self.cfg[idx];
    }
}

struct Observables {
    energy: Vec<f64>,
    magnetisation: Vec<f64>,
}

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

fn magnetisation(cfg: &Configuration) -> f64 {
    return cfg.cfg.iter().sum::<i32>() as f64 / LATSIZE as f64;
}

fn delta_e(cfg: &Configuration, idx: usize) -> i32 {
    return 2*cfg[idx] * (cfg[cfg.neighbours[4*idx]]
                         + cfg[cfg.neighbours[4*idx+1]]
                         + cfg[cfg.neighbours[4*idx+2]]
                         + cfg[cfg.neighbours[4*idx+3]]);
}

fn evolve(cfg: &mut Configuration, energy: &mut f64, beta: f64,
          rng: &mut Rng, nsweep: usize, mut obs: Option<&mut Observables>) -> usize {

    let mut naccept: usize = 0;

    for _sweep in 0..nsweep {
        for _step in 0..LATSIZE {
            let idx = rng.gen_index();

            let delta = delta_e(&cfg, idx);

            if delta <= 0 || (-beta*(delta as f64)).exp() > rng.gen_real() {
                cfg[idx] *= -1;
                *energy += delta as f64;
                naccept += 1;
            }
        }

        if let Some(o) = &mut obs {
            o.energy.push(*energy);
            o.magnetisation.push(magnetisation(&cfg));
        }
    }

    return naccept;
}


fn main() {
    let datadir = Path::new("./data");
    let mut temperatures: Vec<f64> = Vec::new();
    for i in 0..12 {
        temperatures.push((i as f64 + 1.)*0.5);
    }
    prepare_datadir(&datadir, &temperatures);

    let mut rng = Rng::from_seed([138; 32]);

    let mut cfg = Configuration::random(&mut rng);
    let mut energy = hamiltonian(&cfg) as f64;

    let start_time = Instant::now();

    let naccept = evolve(&mut cfg, &mut energy, 1./temperatures[0], &mut rng, NTHERM_INIT, None);
    println!("Initial thermalisation acceptance rate: {}", (naccept as f64)/((NTHERM_INIT*LATSIZE) as f64));

    for (i, temp) in temperatures.iter().enumerate() {
        println!("Running for temperature {}", temp);
        let beta = 1./temp;

        let naccept = evolve(&mut cfg, &mut energy, beta, &mut rng, NTHERM, None);
        println!("  Thermalisation acceptance rate: {}", (naccept as f64)/((NTHERM*LATSIZE) as f64));

        let mut obs = Observables{energy: Vec::new(), magnetisation: Vec::new()};
        let naccept = evolve(&mut cfg, &mut energy, beta, &mut rng, NPROD, Some(&mut obs));
        println!("  Production acceptance rate: {}", naccept as f64 / (NPROD*LATSIZE) as f64);

        write_observables(&datadir.join(format!("{}.dat", i)), &obs);
    }

    let duration = start_time.elapsed();
    println!("Duration in wall clock time: {}", duration.as_secs() as f64
             + (0.001*duration.subsec_millis() as f64));
}
