#!/bin/bash
## Command lines for git

# Help
git config -h # Access help that appears within bash
git config --help # Access help on html

# When using for the first time, set up user name and email in the global setting
git config --global user.name "lingzili"
git config --global user.email "lingzili@bmb.sdu.dk"

# Check your settings
git config --list

# Create a repository
mkdir RNA_Pulldown
cd RNA_Pulldown
git init
ls -a # Check if Git has created a hidden directory called .git

# Steps to track changes by Git
# Step 1: Check the status of the repository
git status
# Step 2: Review changes before adding them to Git
git diff
# Step 3: Add the file to track
git add task_1_fastqc.sh
git status
# Step 4: Commit changes
git commit --message "modify the comments on mkdir for fastqc" # Git will launch the text editor if there is no message

# Change the most recent commit message
git commit --amend -m "New commit message"

# Access change history
git log
git log --all --oneline --graph # Only show git messages

# Access previous commits
git diff HEAD~1 task_1_fastqc.sh

# Access previous commits with unique ID identifier
git diff 9eec9dd606133ac67d125e5cc6dac69f96244ef0 task_1_fastqc.sh
git diff 9eec9d task_1_fastqc.sh

# Access previous commits with the commit message
git show HEAD~1 task_1_fastqc.sh
git show 9eec9d task_1_fastqc.sh

# Acess changes between two versions
git diff 8f0ea0 9eec9d task_1_fastqc.sh

# If you did something wrong and have NOT commit yet, you can replace local changes by git checkout
git checkout 8f0ea0 task_1_fastqc.sh

# Push an existing repository from the command line to GitHub
git remote add origin https://github.com/lingzili/RNA_Pulldown.git
git push -u origin master

# List the remote repository
git remote -v

# Create a local copy of a remote repository (download from GitHub)
git clone ssh://git@github.com/[username]/[repository-name].git	

# List branches (the asterisk denotes the current branch)
git branch	

# Create a new branch
git branch [branch name]

# Swtich to a branch
git checkout [branch name]

# Create a new branch and switch to it
git checkout -b [branch name]

# Merge a branch into a target branch
git merge [source branch] [target branch]

# Delete a branch
git branch -d [branch name]