
import os, string, warnings
from subprocess import *

warnings.simplefilter('once', UserWarning)


repoDir = os.getcwd()

def command(x):
    return str(Popen(x.split(' '), stdout=PIPE).communicate()[0])

def rm_empty(L): return [l for l in L if (l and l!="")]

def getUntracked():
    os.chdir(repoDir)
    status = command("git status")
    if "# Untracked files:" in status:
        untf = status.split("# Untracked files:")[1][1:].split("\n")
        return rm_empty([x[2:] for x in untf if string.strip(x) != "#" and x.startswith("#\t")])
    else:
        return []

def getNew():
    os.chdir(repoDir)
    status = command("git status").split("\n")
    return [x[14:] for x in status if x.startswith("#\tnew file:   ")]

def getModified():
    os.chdir(repoDir)
    status = command("git status").split("\n")
    return [x[14:] for x in status if x.startswith("#\tmodified:   ")]



def gitAdd(fileName, repoDir=repoDir):
    cmd = 'git add ' + fileName
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    print out
    pipe.wait()
    return 

def gitStatus(repoDir=repoDir):
    cmd = 'git status '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    print out,error
    pipe.wait()
    return 

def gitCommit(commitMessage, repoDir=repoDir):
    cmd = 'git commit -m "%s"'%commitMessage
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    print out
    pipe.wait()
    return 

def gitPush(repoDir=repoDir):
    cmd = 'git push '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    pipe.wait()
    return 
    
def gitInit():
    warnings.warn('Thanks to Stack-Overflow users Shane Geiger and Billy Jin for the git wrappers code')
    if os.path.isdir(repoDir + '/.git'):
        warnings.warn('A git repository exists in '+repoDir)
    else:
        cmd = 'git init '
        pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
        (out, error) = pipe.communicate()
        pipe.wait()
        return 
    
def gitLog(repoDir=repoDir):
    cmd = 'git log --all '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    print out,error
    pipe.wait()
    return 