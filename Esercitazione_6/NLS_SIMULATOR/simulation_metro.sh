#!/bin/bash

#METROPOLIS, h=0

mkdir -p "./RESULTS/METROPOLIS"

INPUT_FILE="./INPUT/input.dat"
sed -i "1s/.*/SIMULATION_TYPE        2    1.0     0.0/" "$INPUT_FILE"

cp "./INPUT/CONFIG/config.ising" "./INPUT/CONFIG/config.xyz"

echo "METROPOLIS SIMULATIONS: processing..."

cat > "./INPUT/properties.dat" << EOF
TOTAL_ENERGY
MAGNETIZATION
SPECIFIC_HEAT
SUSCEPTIBILITY

ENDPROPERTIES

/////////////////////////////////////

Possible properties:

LJ Molecular Dynamics or Monte Carlo:

TOTAL_ENERGY
POTENTIAL_ENERGY
KINETIC_ENERGY
TEMPERATURE
PRESSURE
GOFR              100

Only LJ Molecular Dynamics:

POFV              30

Ising 1D MRT^2 or Gibbs:

TOTAL_ENERGY
MAGNETIZATION
SPECIFIC_HEAT
SUSCEPTIBILITY

EOF

echo "------------------ CV U CHI -------------------"

for T in $(seq 2.0 -0.1 0.5)
do
	echo "-----------------------------------------------"
	echo "Running simulation for TEMP=$T"

	sed -i "s/^TEMP.*/TEMP                   $T/" "$INPUT_FILE"
	
	echo "Starting simulation for TEMP=$T"
	
	cd SOURCE
	./simulator.exe

	if [ $? -ne 0 ]; then
		echo "Critical error!"
		exit 1
	fi
	
	cd ..

	echo "Simulation completed, TEMP=$T"
	
	echo "-----------------------------------------------"
	echo "Updating configuration for the next simulation"

	ORIGIN="./OUTPUT/CONFIG/config_20.xyz"
	DESTINATION="./INPUT/CONFIG/config.xyz"

	if [ -f "$ORIGIN" ]; then
		cp "$ORIGIN" "$DESTINATION"
		echo "config.xyz updated successfully"
	else
		echo "Error: $ORIGIN not found"
		exit 1
	fi

	mkdir -p "./RESULTS/METROPOLIS/T_$T"
	mv "./OUTPUT/CONFIG/config_20.xyz" "./RESULTS/METROPOLIS/T_$T/config.xyz"
	mv "./OUTPUT/total_energy.dat" "./RESULTS/METROPOLIS/T_$T/total_energy.dat"
	mv "./OUTPUT/susceptibility.dat" "./RESULTS/METROPOLIS/T_$T/susceptibility.dat"
	mv "./OUTPUT/specific_heat.dat" "./RESULTS/METROPOLIS/T_$T/specific_heat.dat"

	cd SOURCE
	make remove
	cd ..
done

#METROPOLIS, h=0.02

sed -i "1s/.*/SIMULATION_TYPE        2    1.0     0.02/" "$INPUT_FILE"

cp "./INPUT/CONFIG/config.ising" "./INPUT/CONFIG/config.xyz"

cat > "./INPUT/properties.dat" << EOF
MAGNETIZATION

ENDPROPERTIES

/////////////////////////////////////

Possible properties:

LJ Molecular Dynamics or Monte Carlo:

TOTAL_ENERGY
POTENTIAL_ENERGY
KINETIC_ENERGY
TEMPERATURE
PRESSURE
GOFR              100

Only LJ Molecular Dynamics:

POFV              30

Ising 1D MRT^2 or Gibbs:

TOTAL_ENERGY
MAGNETIZATION
SPECIFIC_HEAT
SUSCEPTIBILITY

EOF

echo "---------------- MAGNETIZATION ----------------"

for T in $(seq 2.0 -0.1 0.5) 
do
        echo "-----------------------------------------------"
        echo "Running simulation for TEMP=$T"

        sed -i "s/^TEMP.*/TEMP                   $T/" "$INPUT_FILE"

        echo "Starting simulation for TEMP=$T"

        cd SOURCE
        ./simulator.exe

        if [ $? -ne 0 ]; then
                echo "Critical error!"
                exit 1
        fi

        cd ..

        echo "Simulation completed, TEMP=$T"

        echo "-----------------------------------------------"
        echo "Updating configuration for the next simulation"

        ORIGIN="./OUTPUT/CONFIG/config_20.xyz"
        DESTINATION="./INPUT/CONFIG/config.xyz"

        if [ -f "$ORIGIN" ]; then
                cp "$ORIGIN" "$DESTINATION"
                echo "config.xyz updated successfully"
        else
                echo "Error: $ORIGIN not found"
                exit 1
        fi

        mkdir -p "./RESULTS/METROPOLIS/MAGNET/T_$T"
        mv "./OUTPUT/CONFIG/config_20.xyz" "./RESULTS/METROPOLIS/MAGNET/T_$T/config.xyz"
        mv "./OUTPUT/magnetization.dat" "./RESULTS/METROPOLIS/MAGNET/T_$T/magnetization.dat"

        cd SOURCE
        make remove
        cd ..
done

echo "METROPOLIS SIMULATIONS: completed!"
