# Five dimensional Dirac equation over the reals
## Theoretical background
### Long story made short
Starting point is an article from P Rowlands available on arxiv : http://arxiv.org/abs/quant-ph/0301071

Some articles from José B. Almeida are on Arxiv and give more hints. https://arxiv.org/search/?searchtype=author&query=Almeida%2C+J+B

And articles from Basil Hiley give the link to clifford density element to combine left and right operators. http://www.bbk.ac.uk/tpru/BasilHiley/Bohm-Vienna.pdf

### Details and goals
1. Is it possible to write the Dirac equation with geometric algebra ?
2. Is it possible to write the Dirac equation with a geometric algebra over the reals and NOT the complex ?
3. Is it possible to find solutions of the Dirac equation without the classical physical assumptions (eg. mass is a positive real) ?

Then,
1. Are results from the classic Dirac equation still valid ?
2. Are there any new results with this approach ?
## Technical stack
- The CAS library is sympy,
- the geometric algebra library is galgebra,
- language is python2.
### Run with Jupyter notebook webserver on AWS EC2
This is the preferred method, computations are run on the server, only rendering is happening on the local machine.
An AWS EC2 instance is properly setup with Docker and Git, ssh connection on the instance is open.
```
git clone https://github.com/jdekozak/dirac5d.git
cd dirac5d
docker build --rm -t jupyter/galgebra-notebook .
docker run -p 80:8888 -e NB_UID=500 -e NB_GID=500 --user root -v /home/ec2-user/dirac5d:/home/jovyan/work jupyter/galgebra-notebook
```
### Run with Jupyter notebook webserver on Mac OS X
Computations and rendering are run the local machine.

Docker must be installed https://docs.docker.com/docker-for-mac/install/

Git must be installed https://git-scm.com/download/mac
```
git clone https://github.com/jdekozak/dirac5d.git
cd dirac5d
docker build --rm -t jupyter/galgebra-notebook .
docker run -p 80:8888 -e NB_UID=500 -e NB_GID=500 --user root -v /home/<user>/dirac5d:/home/jovyan/work jupyter/galgebra-notebook
```
Open a browser to your localhost
### Run from a command line terminal
This is the hard way, computations are run locally and all at once, galgebra package generates a PDF document to be viewed with another viewer tool.
```
git clone https://github.com/jdekozak/dirac5d.git
cd dirac5d
git submodule init
git submodule update
python ./dirac5d.py
```

