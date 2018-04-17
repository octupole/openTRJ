'''
Created on Dec 14, 2017

@author: marchi
'''
import sys
import keyring
import getpass
import paramiko
import time

class openFile(object):
    '''
    Open a remote/local file
    '''

    def __init__(self, filename=None,host=None,user=None):
        '''
        Constructor
        '''
        if host:
            if not user:
                print('\nClass "%s": %s \n' % (type(self).__name__,'Username missing in initialization!'))
                sys.exit(1)
            
            password=keyring.get_password(host,user)
            if not password:
                mypass=getpass.getpass()
                keyring.set_password(host,user,mypass)
            password=keyring.get_password(host,user)

            i=1
            while True:
                print("Trying to connect to %s (%i/30)" % (host, i))
                try:
                    client=paramiko.SSHClient()
                    client.load_system_host_keys()
                    client.set_missing_host_key_policy(paramiko.WarningPolicy)
                    client.connect(host,username='marchi',password=password)            
                    print("Connected to %s" % host)
                    break
            
                except paramiko.AuthenticationException:
                    print("Authentication failed when connecting to %s" % host)
                    sys.exit(1)

                except:
                    print("Could not SSH to %s, waiting for it to start" % host)
                    i += 1
                    time.sleep(2)
    # If we could not connect within time limit
                    if i == 30:
                        print("Could not connect to %s. Giving up" % host)
                        sys.exit(1)
            sftp_client = client.open_sftp()
            self.my_file=sftp_client.open(filename)
            print('\nDone remote open\n')
        else:
            self.my_file=open(filename,'r')
            print('\nDone local open\n')

            
    def fp(self):
        return self.my_file

if __name__ == "__main__":

    filename='/ccc/store/cont003/dsv/marchi/AOT-IsoUA/RM10-AML/t600-1100.gyr'
    host='cobalt.ccc.cea.fr'
    myRg=openFile(filename,host,'marchi')
