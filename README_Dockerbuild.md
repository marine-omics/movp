# Docker and Singularity

All of the dependencies for `movp` are provided in a single docker container. 

I the container as follows and then push it to dockerhub. Most of the time this is all that is needed

```bash
docker build -t iracooke/movp .
```


On a systems where `singularity pull` is not working for some readon it might be necessary to manually export as a `.sif` file.  

First on a system with docker
```bash
docker save iracooke/movp -o movp.tar
```
Then copy to a system with singularity installed and then run
```bash
singularity build movp.sif docker-archive://movp.tar
```

