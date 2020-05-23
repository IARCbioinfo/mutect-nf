#!/bin/bash
cd ~/project/
git config --global user.email "alcalan@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://hub.docker.com/api/build/v1/source/33e58fde-38d2-40c6-b7d2-911a0fb2ab5b/trigger/848d00e1-4e32-498a-aaaa-4e17029cdf3c/call/
