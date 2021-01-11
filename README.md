Processing of the datasets for the Blue Cloud sub-contract with VLIZ.

Quick test on BlueCloud

```bash
rm -R bluecloud-plankton-master.tar.gz bluecloud-plankton-master
wget https://gitlab.uliege.be/gher/bluecloud-plankton/-/archive/master/bluecloud-plankton-master.tar.gz
tar -xzf bluecloud-plankton-master.tar.gz
cd bluecloud-plankton-master
julia --project=@. -e "import Pkg; Pkg.test()
```