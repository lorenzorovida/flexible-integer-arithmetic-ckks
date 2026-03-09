# A flexible and polynomial framework for integer arithmetic in CKKS

[![CMake on multiple platforms](https://github.com/lorenzorovida/flexible-integer-arithmetic-ckks/actions/workflows/cmake-multiple-platform.yml/badge.svg)](https://github.com/lorenzorovida/flexible-integer-arithmetic-ckks/actions/workflows/cmake-multiple-platform.yml)
<br>
<a href="https://eprint.iacr.org/2025/1150"><img src="imgs/preprint_icon.svg" alt="Link to the preprint PDF" ></a>

<img src="imgs/console.gif" width="650"/>

---

This repository contains an OpenFHE-based project that implements a computational device that handles power-of-two exact arithmetic, up to 256-bits.

The reference paper for this work is [A flexible and polynomial framework for integer arithmetic in CKKS](https://ia.cr/2026/450).

| :exclamation: Keep in mind!                                  |
|:-------------------------------------------------------------|
| This repo is still WIP, so feel free to report any issue! ^^ |

## Fresh start
One of the key ideas behind this work is *simplicity*. It can indeed be used by simply installing a custom OpenFHE fork and by compiling the CMake project.

We require a custom fork of OpenFHE (that is updated to v1.5.0 and it includes the StC-first bootstrapping), this allows to use the two functions `EvalChebyshevSeriesPSBatchRepeated` and `EvalBootstrapStCFirstBits`. The first generalizes the functionality to evaluate a Chebyshev polynomial over a ciphertext to evaluate $n$ Chebyshev polynomials, one for each slot of the ciphertext. The second allows to evaluate a cleaning bootstrapping operation à la [[BCKS24]](https://eprint.iacr.org/2024/767).


1) Install the `repeated_poly_and_stcboot` branch from [this](https://github.com/lorenzorovida/openfhe-development-chebyshevSIMD) custom fork of OpenFHE 

```
git clone --branch repeated_poly_and_stcboot --single-branch https://github.com/lorenzorovida/openfhe-development-chebyshevSIMD
cd openfhe-development-chebyshevSIMD
mkdir build
cmake ..
sudo make install
```


2. Now, navigate to your desired folder and install this repository

```
git clone https://github.com/lorenzorovida/flexible-integer-arithmetic-ckks
cd flexible-integer-arithmetic-ckks
mkdir build
```

The `CMakeLists.txt` **must** point to the custom version of OpenFHE. **Be sure** to replace all the instances of `/Users/myuser/openfhe-development-chebyshevSIMD/build` in the `CMakeLists.txt` file to your custom build folder! After that, build the repository:
```
cmake ..
make
```



> [!IMPORTANT]  
> The compilation will fail if you are using the original version of OpenFHE, use the fork!


3. Finally, you can test if everything works by running

```
./FlexibleIntsCKKS --test
```


## Custom parameters
There are three parameters that can be passed to `./FlexibleIntsCKKS`

1) `--ring <size>`

Sets the (logarithm of the) ring size. `<size>` must be an integer in (12, 13, 14, 15, 16).
Example:
```
./FlexibleIntsCKKS --ring 16
```

2) `--bits <bits>` 

Sets the word size (number of bits per word). Only the following values are supported: 8, 16, 32, 64, 128, 256.
If an unsupported value is provided, the program will display an error.
Example:
```
./FlexibleIntsCKKS --bits 64
```

3) `--verbose <value>` 

Sets the verbosity level of program output. `<value>` must be either 0, 1 or 2. Higher numbers mean more detailed output.
Example:
```
./FlexibleIntsCKKS --verbose 2
```

4) `--input`

Sets the program in input mode, i.e., the program will ask you to set the number of desired bits and to manually insert the values you want to compute operations on.

## Examples of usage

#### Toy parameters
You can test 64-bits operation, in the $N=2^{12}$ ring, with the following code:
```
./FlexibleIntsCKKS --ring 12 --bits 64
```
This allows to test the code in a unsecure environment, useful if playing around with the code.

#### More secure parameters
You can test 64-bits operation, in the $N=2^{16}$ ring, with the following code:
```
./FlexibleIntsCKKS --ring 16 --bits 64
```

#### Input mode
You can use the program and manually insert some values to operate on
```
./FlexibleIntsCKKS --ring 12 --input
```
