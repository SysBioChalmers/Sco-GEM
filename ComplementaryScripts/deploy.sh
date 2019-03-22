#!/usr/bin/env bash

# do NOT set -v or your GitHub API token will be leaked!
set -e # exit with nonzero exit code if anything fails

source_branch="chore/travis"

if [[ "${TRAVIS_PULL_REQUEST}" != "false" || "${TRAVIS_BRANCH}" != "${source_branch}" || "${TRAVIS_REPO_SLUG}" != "SysBioChalmers/Sco-GEM" ]]; then
    echo "Skip deploy."
    exit 0
else
  echo "Starting deploy to ${DEPLOY_BRANCH}..."
fi

# configure git
git config --global user.email "deploy@travis-ci.org"
git config --global user.name "Travis CI Deployment Bot"

git status
git checkout -- ./ComplementaryScripts/deploy.sh
memote run --exclusive test_basic #--skip-unchanged
git status
git checkout ${DEPLOY_BRANCH}
ls
memote --directory="./results/" --filename="index.html"
git status
git add "index.html"
git commit -m "Travis build ${TRAVIS_BUILD_NUMBER}"
git push




# # clone the deploy branch
# cd "${HOME}"
# git clone --quiet --branch=${DEPLOY_BRANCH} https://${GITHUB_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git ${DEPLOY_BRANCH} > /dev/null

# # copy the results from the current memote run to deploy dir
# cp "${TRAVIS_BUILD_DIR}/results/${TRAVIS_COMMIT}.json.gz" "${HOME}/${DEPLOY_BRANCH}/results/"

# # create the report pointing to the history stored in deploy branch
# # need to be in build directory to access git history
# cd "${TRAVIS_BUILD_DIR}"
# memote --directory="${HOME}/${DEPLOY_BRANCH}/results/" --filename="${HOME}/${DEPLOY_BRANCH}/index.html" report

# #add, commit and push files
# cd "${HOME}/${DEPLOY_BRANCH}"
# git add "results/${TRAVIS_COMMIT}.json.gz"
# git add "index.html"
# git commit -m "Travis build ${TRAVIS_BUILD_NUMBER}"
# git push --quiet origin "${DEPLOY_BRANCH}" > /dev/null

echo "Done."
