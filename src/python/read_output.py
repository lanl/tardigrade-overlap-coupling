import os

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xml.parsers.expat

def get_attribute_data(xml_node, path='.'):
    
    if ("Attribute" != xml_node.tag):
        raise ValueError("XML node is not of type 'Attribute'")
        
    value = None
        
    for child in xml_node:
                            
        if "DataItem" == child.tag:
                                
            value = read_data(child, path=path)
                
            break
            
    return value

def get_set(xml_node, path='.'):
    
    if ("Set" != xml_node.tag):
        raise ValueError("XML node is not of type 'Set'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break     
        
    return value

def get_geometry(xml_node, path='.'):
    
    if ("Geometry" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Geometry'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break
            
    return value

def get_topology(xml_node, path='.'):
    
    if ("Topology" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Topology'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break
            
    return value

def read_data(xml_node, path='.'):
    
    if ("DataItem" != xml_node.tag):
        raise ValueError("XML node is not of type 'DataItem'")
        
    shape = [int(v) for v in xml_node.attrib["Dimensions"].split()]
    
    if xml_node.attrib["DataType"] == "Float":
        dtype = float
    else:
        dtype = int
        
    if xml_node.attrib["Format"] == "XML":
        value = np.hstack([np.array(line.split()).astype(dtype) for line in xml_node.text.split("\n")]).reshape(shape)
        
    else:
        
        filename, keyname = xml_node.text.split(':')
        
        with h5py.File(os.path.join(path, filename), 'r') as h5file:
            
            value = np.array(h5file[keyname]).astype(dtype).reshape(shape)
        
    return value