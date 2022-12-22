# pb-hifi-viral

WDL workflows and dockerfiles for tools related to PacBio's HiFi Viral SARS-CoV-2 sequence analysis

This repo is under active development and will change without notice.

## Docker

This repo includes a Dockerfile, which is used for building the docker image that will be used by the `run_sample` WDL task in the WDL workflow `pb_sars_cov2_kit.wdl`

To build the docker image locally, run the following:

```bash
# clone the code locally
git clone https://github.com/theiagen/pb-hifi-viral.git
cd pb-hifi-viral

# build the docker image on your local machine, build both the app and test layer
docker build --target test -t theiagen/pb-hifi-viral:latest docker/

# build the docker image on your local machine, only build app layer, NOT test layer
# Use this for deploying new docker image to quay
docker build --target app -t quay.io/theiagen/pb-hifi-viral:latest docker/

# (optionally) push the docker image to Theiagen quay repo for usage in Terra
docker push quay.io/theiagen/pb-hifi-viral:latest
```
