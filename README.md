## Ratlas app

__Purpose of the repository__: 

The goal of this repository is to serve as a central location to all who will be working in the app to avoid effort duplication and to follow best practices for version control.

For this to work constructively: 

1) `git clone` this repository to your local computer and feel free to make changes to the repository. If you would like to contribute to this repository please follow the following git flow: https://levelup.gitconnected.com/semantic-versioning-with-git-flow-and-the-marvelous-way-to-go-there-b9f97b90455c?gi=3901bb647b71 ; for more information on this please contact a member from U-BDS. In short, the `main` branch is the stable version of the app, and new features should be committed to the `development` branch (ideally by first creating `feature_*` branch from `development` first). Once the newly developed feature is fully tested it may be merged to `main`. U-BDS provides a detailed guide on this protocol (request access if needed).

2) Make sure you have the datasets needed by the app. Download them from https://uab.box.com/s/jc4bgyqmzkrol74gtzo2hgcfdkxhikxg (if you need access, ask Lara), and place them in the `lean_datasets` folder from this repository. __You do not need to download the `other_files` sub-directory (these are just files which are present for solidarity or development needs)__.

## Current state:

The git hosting location has now been migrated to GitHub as a brand new repo. To avoid security concerns which may be linked to old commits which date back to September 30th, 2019 (as we are moving forward with plans to host the app in cloud.rc), git mirroring was not performed, and all code in GitHub was migrated as a fresh repository. Thus the first commit in GitHub is representative of all work performed since the original launch up until June 2022. This included, original datasets in original launch (adult NAc, primary striatal culture dataset), addition of VTA dataset, addition of rn7 mode for all datasets aligned to the rn6 reference, improvement in performance or usage of the app (e.g.: searching gene names which do not have a gene name), and initial development work for ATAC-seq datasets.

## Running the app with Docker

A Docker container has also been developed for this app which will eventually be launched in our own cloud.rc rather than shinyapps.io.

In the meantime, you can launch a containarized version of the app by downloading the image which contains the shiny/app dependencies from Docker Hub: 

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

Above directions also apply to launching the app in cloud.rc but note the use of the `--restart always` flag instead of `--rm` due to automatic server reboots:

```bash
sudo docker run -d --restart always --user shiny -p 3838:3838 -v `pwd`/ratlas:/srv/shiny-server/ -v `pwd`/shiny_app_logs:/var/log/shiny-server uabbds/ratlas:latest
```