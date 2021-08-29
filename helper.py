import numpy as np

def constant_refractive_index(n, ws):
    '''
    Args:
        n (float): constant refractive index
        ws (float list): wavenumber
    Return:
        refractive_index (complex array): refractive index
    '''
    refractive_index = complex(n, 0)
    refractive_index = refractive_index + np.array(ws) * 0
    return refractive_index

def set_layer(refractive_index, thickness=0, coherence=True):
    '''
    Args:
        refractive_index (complex array): complex refractive index
        thickness (float): thickness of layer (unit:cm)
        coherence (boolean): coherence in layer
    Return:
        layer information (dictionary)
    '''  
    return {'refractive_index': refractive_index, 'thickness': thickness, 'coherence': coherence}