"""
Benchmark different implementations.

Automatically adjusts lattice size and temperature but does not change
the number of sweeps. They need to be set in the base implementations manually.
"""

from pathlib import Path
import subprocess
import tempfile
import re
import pickle
import shutil

# output file for timings
BENCHMARK_FILE = Path(__file__).parent/"ising.ben"

# source files
PYTHON_SRC = Path(__file__).parent/"../python/ising.py"
CPP_SRC = Path(__file__).parent/"../cpp/base/ising.cpp"
RUST_SRC = Path(__file__).parent/"../rust/src/main.rs"


def run_python(sizes, wdir):
    "Run the Python implementation."

    # regexes to extract and replace stuff
    nxre = re.compile(r"NX = (\d+).*$")
    nyre = re.compile(r"NY = (\d+).*$")
    tempre = re.compile(r"TEMPERATURES = ([^\#]+)(\#.*)?$")
    timere = re.compile(r"Duration in wall clock time: ([\d\.]+)s")

    # read soruce code
    with open(PYTHON_SRC, "r") as srcf:
        src = srcf.read()

    # script to run
    script = wdir/"ising.py"

    times = []
    for size in sizes:
        # insert size and fix temperature
        src = "\n".join(re.sub(nxre, f"NX = {size}",
                               re.sub(nyre, f"NY = {size}",
                                      re.sub(tempre, "TEMPERATURES = [1.]", line)))
                        for line in src.split("\n"))
        with open(script, "w") as outf:
            outf.write(src)

        # run the benchmark
        proc = subprocess.run(["python3", str(script), str(wdir/"data")],
                              capture_output=True)

        # extract run time
        for line in proc.stdout.decode("utf-8").split("\n"):
            match = re.match(timere, line)
            if match:
                times.append(float(match[1]))
                break

    return times

def run_cpp(sizes, wdir):
    "Run the C++ implementation."

    # regexes to extract and replace stuff
    nxre = re.compile(r"constexpr Index NX = (\d+);.*$")
    nyre = re.compile(r"constexpr Index NY = (\d+);.*$")
    timere = re.compile(r"Duration in wall clock time: ([\d\.]+)s")

    # read source code
    with open(CPP_SRC, "r") as srcf:
        src = srcf.read()

    # temporary files to run the program
    cmake_file_in = CPP_SRC.parent/"CMakeLists.txt"
    cmake_file = wdir/"CMakeLists.txt"
    build_dir = wdir/"build"
    source = wdir/"ising.cpp"

    # copy over cmake file and run cmake
    shutil.copy(cmake_file_in, cmake_file)
    shutil.copy(CPP_SRC, source)
    build_dir.mkdir()
    subprocess.run(["cmake", "..", "-DCMAKE_BUILD_TYPE=RELEASE"], cwd=build_dir)

    times = []
    for size in sizes:
        # insert size and fix temperature
        src = "\n".join(re.sub(nxre, f"constexpr Index NX = {size};",
                               re.sub(nyre, f"constexpr Index NY = {size};",
                                      line.replace("constexpr auto temperatures = listTemperatures();",
                                                   "constexpr std::array<double,1> temperatures{1.,};")))
                        for line in src.split("\n"))
        # write source and compile
        with open(source, "w") as outf:
            outf.write(src)
        subprocess.run(["make"], cwd=build_dir)

        # run the program
        proc = subprocess.run(["ising", str(wdir/"data")], cwd=build_dir,
                              capture_output=True)

        # extract run time
        for line in proc.stdout.decode("utf-8").split("\n"):
            match = re.match(timere, line)
            if match:
                times.append(float(match[1]))
                break
    return times

def run_rust(sizes, wdir):
    "Run the Rust implementation."

    # regexes to extract and replace stuff
    nxre = re.compile(r"const NX: usize = (\d+);.*$")
    nyre = re.compile(r"const NY: usize = (\d+);.*$")
    timere = re.compile(r"Duration in wall clock time: ([\d\.]+)s")

    # read source code
    with open(RUST_SRC, "r") as srcf:
        src = srcf.read()

    # make backup of source file
    backup = wdir/"main.rs"
    shutil.copy(RUST_SRC, backup)

    times = []
    for size in sizes:
        # insert size and fix temperature
        src = "\n".join(re.sub(nxre, f"const NX: usize = {size};",
                               re.sub(nyre, f"const NY: usize = {size};",
                                      line.replace("let temperatures = list_temperatures();",
                                                   "let temperatures = vec!(1.);")))
                        for line in src.split("\n"))
        # write source and compile
        with open(RUST_SRC, "w") as outf:
            outf.write(src)

        # run the program
        proc = subprocess.run(["cargo", "run", str(wdir/"data")], cwd=RUST_SRC.parent.parent,
                              capture_output=True)

        # extract run time
        for line in proc.stdout.decode("utf-8").split("\n"):
            match = re.match(timere, line)
            if match:
                times.append(float(match[1]))
                break

    # restore original source
    shutil.copy(backup, RUST_SRC)

    return times


def main():
    "Run the benchmarks."

    sizes = [4, 16, 24, 32]

    times = dict()
    with tempfile.TemporaryDirectory() as wdir:
        wdir = Path(wdir)
        times["python"] = run_python(sizes, wdir)
        times["c++"] = run_cpp(sizes, wdir)
        times["rust"] = run_rust(sizes, wdir)

    # save benchmark to file
    pickle.dump({"xlabel": "Nx*Ny",
                 "ylabel": "time / s",
                 "xvalues": [size**2 for size in sizes],
                 "yvalues": times},
                open(BENCHMARK_FILE, "wb"))



if __name__ == "__main__":
    main()
