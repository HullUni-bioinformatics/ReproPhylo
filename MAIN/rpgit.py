
import os, string, warnings, time

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



def gitAdd(fileName,  repoDir=repoDir):
    cmd = 'git add ' + fileName
    #print 'DEBUG gitAdd: %s'%cmd
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    #print 'DEBUG gitAdd STDOUT: %s'%out
    #print 'DEBUG gitAdd STDERR: %s'%error
    pipe.wait()
    return out, error

def gitStatus(repoDir=repoDir):
    cmd = 'git status '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    pipe.wait()
    return out, error

def gitCommit(commitMessage, repoDir=repoDir):
    cmd = 'git commit -m "%s"'%(commitMessage.replace('\"','\''))
    #print 'DEBUG gitCommit: %s'%cmd
    try:
        pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
        (out, error) = pipe.communicate()
    except:
        try:
            hndl = open('faulty_git_commit_messages', 'a')
            hndl.write(commitMessage + '\n')
            hndl.close()
            commitMessage = commitMessage.splitlintes()[0] + '\nSee faulty_git_commit_messages for more details.'
            cmd = 'git commit -m "%s"'%(commitMessage.replace('\"','\''))
            pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
            (out, error) = pipe.communicate()
        except:
            hndl = open('faulty_git_commit_messages', 'a')
            hndl.write(commitMessage + '\n')
            hndl.close()
            commitMessage = "Commit successful but something very wrong with message.\nSee faulty_git_commit_messages for more details."
        cmd = 'git commit -m "%s"'%(commitMessage.replace('\"','\''))
        pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
        (out, error) = pipe.communicate()
        
    #print 'DEBUG gitCommit STDOUT: %s'%out
    #print 'DEBUG gitCommit STDERR: %s'%error
    pipe.wait()
    return out, error

def gitPush(repoDir=repoDir):
    cmd = 'git push '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    pipe.wait()
    return out, error
    
def gitInit():
    
    warnings.warn('Thanks to Stack-Overflow users Shane Geiger and Billy Jin for the git wrappers code')
    if os.path.isdir(repoDir + '/.git'):
        warnings.warn('A git repository exists in '+repoDir)
        return ('A git repository exists in '+repoDir, 'None')
    else:
        emailcmd = 'git config --global user.email'
        pipe = Popen(emailcmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
        (emailout, emailerror) = pipe.communicate()
        if not emailout or (emailout and not '@' in emailout):
            raise RuntimeError('Git: set your email with \'!git config --global user.email "your_email@example.com"\''+
                              ' or disable git (the ! is needed in IPython Notebook. In a terminal, ommit it)')
        pipe.wait()
        cmd = 'git init '
        pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
        (out, error) = pipe.communicate()
        pipe.wait()
        warnings.warn('A git repository was created in %s.'%repoDir)
        return out, error
    
def gitLog(repoDir=repoDir):
    cmd = 'git log --all '
    pipe = Popen(cmd, shell=True, cwd=repoDir,stdout = PIPE,stderr = PIPE )
    (out, error) = pipe.communicate()
    pipe.wait()
    return out,error