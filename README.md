# CH502:: Quantum Chemistry

## Hartree-Fock Calculation of H2 Molecule

### 0. Environmental Settings
```
    pip install -r requirements.txt
```

### 1. Run HF with .xyz file
```
    python main.py $PATH_TO_H2.xyz
```

### 2. Run HF to draw Dissociation Curve
```
    python plot.py > /dev/null
```

<img src="https://github.com/WonhoZhung/CH502/blob/main/Hartree-Fock_H2/Dissociation_Curve_H2.png" height="400">
<img src="https://github.com/WonhoZhung/CH502/blob/main/Hartree-Fock_H2/Energy_Curve_H2.png" height="400">


### Credits ↓


                             ..:::::::::::::::....
                   .::== .#%*  +++++++++++++++++++: =##= .:..
           .::=+++++++: =@@@  +++++++++++++++++++. #@@# .+++++==:.
         .:=+++++++++++:    .=++++++++++++++++++++.    .+++++++++=:
               ..:==++++++++++++++++++++++++++++++++++++++=::.
                       ..:::====++++++++++++++===:::..
                             .:+**************
                         .:++++=. +++++++++=+++:
                  .+++=++++=.     +++++++++= :+++.
                  +++++:.         +++++++++=   =+++:=:
                   .::.           +++++++++=     :+++++
                                  ++=::::++=      :+++:
                                  ++:    ++=
                                  ++:    ++=
                             .::::++:    ++=.:::
                            .*++++++:    +++++++*
                            
