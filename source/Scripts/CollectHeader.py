import os
import string

file_name = 'lib.h'
folder='c:/Projects/bFluid/source/lib'

if os.path.exists(file_name):
    os.remove(file_name)

fo=file(file_name,'w')

folder='c:/Projects/bFluid/source/lib'

for roots,dirs,files in os.walk(folder):
    for header in files:
        if header.endswith(('.h','.hpp')):
            header=str(os.path.join(roots,header))
            header_new = header.lstrip(folder)
            x=header_new.replace('\\','/')
            y=str('#include <'+x+'>\n')
            y=y.replace('</','<')
            fo.write(y)
fo.close()

print 'lib.h is generated!'


            
