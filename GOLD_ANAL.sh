WORKDIR=/home/mspencer/Documents/Corrections/Au/Kinematics


for k in $(seq 37 43);
do
number=$(( k*25 ))
./kinan "$number"_kin.dat $number
done
