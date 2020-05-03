
if [ $# -eq 0 ]
  then
    echo "No name supplied"
else
  echo $1

  conda create -n $1 python=3. -y
  conda install -n $1 numpy -y
  conda install -n $1 rdkit -c rdkit -y
  conda install -n $1 openbabel -c openbabel -y
  conda install -n $1 pytorch -c pytorch -y
  conda install -n $1 scipy -y
  conda install -n $1 pytest -y
  conda install -n $1 tqdm -y
  conda install -n $1 -c conda-forge bayesian-optimization
  conda install -n $1 -c conda-forge ase
  source activate $1
  pip install qml

fi
