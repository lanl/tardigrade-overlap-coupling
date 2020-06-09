import yaml

document = """
coupling:
    micromorphic:
        filename: micromorphic_output_file.xdmf
        ghost blocks: 0
        free blocks: 1
    dns:
        filename: dns_output_file.xdmf
        ghost blocks: 0
        free blocks: 1
    linkage:
        micromorphic ghost to dns ghost:
        micromorphic free to dns ghost:
        micromorphic ghost to dns free:
        micromorphic cree to dns free: 
"""

print( yaml.full_load( document ) )
