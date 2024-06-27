# Ratlas app

## Data location (for development purposes)

This section is written for the app developer(s) or contributors. **General users are encouraged to access the web app by the public instance / URL: <https://ratlas.org/>**

**The Ratlas datasets** The on-disk data can be downloaded from the following link <https://uab.box.com/s/jc4bgyqmzkrol74gtzo2hgcfdkxhikxg>. Only the lab and any developers have access to this folder. The folder has a `README` with relevant information about the data versioning.

If you would like access to the original Seurat object, please contact Dr. Jeremy Day. You may find his contact information at the lab web page: <https://day-lab.org/uab>

## Running the app with Docker

A Docker container has been created for the deployment of the app in UAB's cloud.rc.

To make use of the container in your own computer (for local tests etc.), please see the following steps:

* `git clone` this app repository and `cd` to it from your computer

* Run the container with the following command (where `pwd` should be path upstream of the ratlas directory where app code is located. Thus, be sure to go to this cloned repository):

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

The above directions also apply to launching the app in cloud.rc. But please note that you should deploy the app with the use of the `--restart always` flag instead of `--rm` due to automatic server reboots:

```bash
sudo docker run -d --restart always --user shiny -p 3838:3838 -v `pwd`/ratlas:/srv/shiny-server/ -v `pwd`/shiny_app_logs:/var/log/shiny-server uabbds/ratlas:latest
```