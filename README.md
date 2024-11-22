# DpDse
Dynamic state estimator for electrical grids using dynamic phasor formulation


## Install

### Docker
First, you need to install [Docker]([https://link-url-here.org](https://docs.docker.com/get-docker/)). Then, you could build the docker image by yourself as follows:

1. Clone the repository <br>
```
$ git clone git@github.com:your_git_username/DSE.git
```

3. In the repository, there is a Docker file with all required dependencies<br>
```
cd DSE
$ sudo docker build -t dse .
```

3. Next, run a Docker container<br>
```
$ docker run -it -p 8888:8888 dse bash
```

4. Open jupyter lab
```
jupyter lab --ip="0.0.0.0" --allow-root --no-browser
```
