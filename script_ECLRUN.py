#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/python

import subprocess
import os
import sys
import glob
import os.path

ECL_path='./AA222_EclipseFile' #Folder name that contains the Eclipse files, currently on desktop
os.chdir(ECL_path)

# Run Eclipse
#runEclipse = subprocess.Popen(["eclrun","eclipse", ECL_name], 
#                    stdin =subprocess.PIPE,
#                    stdout=subprocess.PIPE,
#                    stderr=subprocess.PIPE,
#                    universal_newlines=True,
#                    bufsize=0)
os.system('cmd /k "eclrun eclipse AA222.data"')