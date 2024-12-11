import numpy
import matplotlib.pyplot as plt 
import json
from jinja2 import Environment, FileSystemLoader




def Create_json_wrapping(ka,kb,r,inter_str):

    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping.txt')
    output_from_parsed_template = template.render(KA = ka, KB = kb,radius = r,xpos = r*1.9 ,interaction=inter_str)

    data = json.loads(output_from_parsed_template)


    with open('../Config_files/Wrapping_strg_{}_radius_{}_KA_{}_KB_{}.json'.format(inter_str,r,ka,kb), 'w') as file:
        json.dump(data, file, indent=4)



Create_json_wrapping(0.05,10.0,0.3,200)

