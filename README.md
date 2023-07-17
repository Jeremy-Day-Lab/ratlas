# Ratlas app

## Data location (for development purposes)

This section is written for the app developer(s) or contributors. **General users are encouraged to access the web app by the public instance / URL: <https://day-lab.shinyapps.io/ratlas/>**

**The Ratlas datasets** The data can be downloaded from the following link <https://uab.box.com/s/jc4bgyqmzkrol74gtzo2hgcfdkxhikxg>. Only the lab and any developers have access to this folder. The folder has a `README` with relevant information about the data versioning.

## Running the app with Docker

A Docker container has also been developed for this app which will eventually be launched in our own cloud.rc rather than shinyapps.io.

In the meantime, for development, you can launch a containarized version of the app by downloading the image which contains the shiny/app dependencies from Docker Hub: 

* `git clone` this app repository and `cd` to it from your computer

* Run the container with the following command (again where `pwd` should be path upstream of the ratlas directory where app code is located. Thus, be sure to go to this cloned repository):

```bash
docker run -d --rm --user shiny -p 3838:3838 -v `pwd`/ratlas:/srv/shiny-server/ -v `pwd`/shiny_app_logs:/var/log/shiny-server uabbds/ratlas:latest
```

Open your browser, and go to local `localhost:3838`

Once finished don't forget to stop the container with

```
docker ps # find docker container id
docker stop <container_id>
```

## Running the app with Docker in cloud.rc

Above directions also apply to launching the app in cloud.rc for testing purposes. Note the use of the `--restart always` flag instead of `--rm` due to automatic server reboots:

```bash
sudo docker run -d --restart always --user shiny -p 3838:3838 -v `pwd`/ratlas:/srv/shiny-server/ -v `pwd`/shiny_app_logs:/var/log/shiny-server uabbds/ratlas:latest
```