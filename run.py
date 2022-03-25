import subprocess as sp
import os
import json

Table =  json.loads(open('Table.json').read())
print(Table['elements'][0]['symbol'])
for elem in Table['elements']:
    print(elem['symbol'])
