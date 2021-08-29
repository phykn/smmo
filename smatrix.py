import numpy as np

def get_cos_angle(n_0, n_i, q_0):
    assert len(n_0) == len(n_i), "size mismatch"
    assert q_0 >= 0 and q_0 < 90, "Please set 0 <= q_0 < 90"
    '''
    Args:
        n_0 (complex array): refractive index of 0th layer
        n_i (complex array): refractive index of ith layer
        q_0 (float): incidence angle of 0th layer
    Return:
        cos_qi (complex array): cos angle of ith layer
    '''
    
    q_0 = np.full(len(n_0), complex(q_0*np.pi/180, 0), dtype=complex)  
    cos_qi = np.sqrt(1 - np.square((n_0 * np.sin(q_0) / n_i)))
    return cos_qi

def get_kz(n_i, cos_qi, ws):
    assert len(n_i) == len(cos_qi) == len(ws), "size mismatch"
    '''
    Args:
        n_i (complex array): refractive index of ith layer
        cos_qi (complex array): cos angle of ith layer
        ws (float array): list of wavenumber
    Return:
        kz (complex array): wavevector in z direction
    '''
    
    kz = 2 * np.pi * n_i * ws * cos_qi
    return kz

def get_fresnel_coeff_ij(n_i, n_j, cos_qi, cos_qj, rt='r', sp='s'):
    assert rt in ['r', 't'], "rt not in ['r', 't']"
    assert sp in ['s', 'p'], "sp not in ['s', 'p']"
    assert len(n_i) == len(n_j) == len(cos_qi) == len(cos_qj), "size mismatch"
    '''
    Args:
        n_i (complex array): refractive index of ith layer
        n_j (complex array): refractive index of jth layer
        cos_qi (complex array): cos angle of ith layer
        cos_qj (complex array): cos angle of jth layer
        rt: 'r' (refelction) or 't' (transmission)
        sp: 's' (s-polarization) or 'p' (p-polarization) 
    Return:
        Fresnel coefficient (complex array): rs, rp, ts, tp
    '''    
    
    if (rt == 'r') and (sp == 's'):   
        rs = (n_i*cos_qi-n_j*cos_qj)/(n_i*cos_qi+n_j*cos_qj)
        return rs
    elif (rt == 'r') and (sp == 'p'):
        rp = (n_j*cos_qi-n_i*cos_qj)/(n_j*cos_qi+n_i*cos_qj)
        return rp
    elif (rt == 't') and (sp == 's'):
        ts = 2*n_i*cos_qi/(n_i*cos_qi+n_j*cos_qj)
        return ts
    elif (rt == 't') and (sp == 'p'):
        tp = 2*n_i*cos_qi/(n_j*cos_qi+n_i*cos_qj)
        return tp
    
def get_s_matrix(layers, n_0, q_0, ws, sp='s'):
    assert sp in ['s', 'p'], "sp not in ['s', 'p']"
    assert len(n_0) == len(ws), "size mismatch"
    '''
    Args:
        layers (dictionary list): list of layers
        n_0 (complex array): refractive index of 0th layer
        q_0 (float): incidence angle of 0th layer
        ws (float array): list of wavenumber
        sp: 's' (s-polarization) or 'p' (p-polarization) 
    Return:
        S_11, S_12, S_21, S_22 (complex array)
    '''    
    
    # Add a pseudo final layer
    layer_f = layers[-1].copy()
    layer_f['thickness'] = 0.0
    layer_f['coherence'] = True
    layers.append(layer_f)

    # Get Numbers
    n_layer = len(layers)
    n_ws    = len(ws)

    # Initialize
    n     = np.zeros((n_layer, n_ws), dtype=complex)
    d     = np.zeros((n_layer)      , dtype=complex)   
    cos_q = np.zeros((n_layer, n_ws), dtype=complex)
    phase = np.zeros((n_layer, n_ws), dtype=complex)    
    for i, layer in enumerate(layers):
        n[i]     = layer['refractive_index']
        d[i]     = layer['thickness']      
        cos_q[i] = get_cos_angle(n_0, n[i], q_0)
        phase[i] = np.exp(1j * get_kz(n[i], cos_q[i], ws) * d[i])
        
    # Get S-matrix
    S_11 = np.full((n_ws), 1, dtype=complex)
    S_12 = np.full((n_ws), 0, dtype=complex)
    S_21 = np.full((n_ws), 0, dtype=complex)
    S_22 = np.full((n_ws), 1, dtype=complex)
    for i in range(n_layer-1):
        j = i + 1        
        t = get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], rt='t', sp=sp)
        r = get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], rt='r', sp=sp)

        I_11 = 1 / t
        I_12 = r / t
        I_21 = I_12
        I_22 = I_11

        S_11 = phase[i] * S_11 / (I_11 - phase[i] * S_12 * I_21)
        S_12 = ((phase[i] * S_12 * I_22 - I_12) * phase[j]) / (I_11 - phase[i] * S_12 * I_21)
        S_21 = S_22 * I_21 * S_11 + S_21
        S_22 = S_22 * I_21 * S_12 + S_22 * I_22 * phase[j]

    return S_11, S_12, S_21, S_22

def get_tr_matrix(layers, n_0, q_0, ws, sp='s'):
    assert sp in ['s', 'p'], "sp not in ['s', 'p']"
    assert len(n_0) == len(ws), "size mismatch"
    '''
    Args:
        layers (dictionary list): list of layers
        n_0 (complex array): refractive index of 0th layer
        q_0 (float): incidence angle of 0th layer
        ws (float array): list of wavenumber
        sp: 's' (s-polarization) or 'p' (p-polarization) 
    Return:
        T_12, R_12, T_21, R_21 (float array)
    '''
    
    forward_layers = layers   
    S_11, S_12, S_21, S_22 = get_s_matrix(forward_layers, n_0, q_0, ws, sp=sp)
    T_12 = (S_11 * np.conjugate(S_11)).real
    R_12 = (S_21 * np.conjugate(S_21)).real

    backward_layers = layers[::-1]
    S_11, S_12, S_21, S_22 = get_s_matrix(backward_layers, n_0, q_0, ws, sp=sp)
    T_21 = (S_11 * np.conjugate(S_11)).real
    R_21 = (S_21 * np.conjugate(S_21)).real

    return T_12, R_12, T_21, R_21

def get_TR(layers, n_0, q_0, ws, sp='s'):    
    assert sp in ['s', 'p'], "sp not in ['s', 'p']"
    assert len(n_0) == len(ws), "size mismatch"
    '''
    Args:
        layers (dictionary list): list of layers
        n_0 (complex array): refractive index of 0th layer
        q_0 (float): incidence angle of 0th layer
        ws (float array): list of wavenumber
        sp: 's' (s-polarization) or 'p' (p-polarization) 
    Return:
        T_total, R_total (float array)
    '''    
    
    pseudo_layers = []
    for i, layer in enumerate(layers):
        if layer['coherence'] == False:
            pseudo_layer = {'refractive_index': layer['refractive_index'], 'thickness': 0.0, 'coherence': True}
            pseudo_layers.append(pseudo_layer)
            pseudo_layers.append(layer)
            pseudo_layers.append(pseudo_layer)
        else:
            pseudo_layers.append(layer)

    blocks = []
    block = pseudo_layers[0:1]
    for i in range(len(pseudo_layers)):
        if i != len(pseudo_layers)-1:
            layer_i = pseudo_layers[i]
            layer_j = pseudo_layers[i+1]  
            if layer_i['coherence'] == layer_j['coherence']:
                block.append(layer_j)
            else:
                blocks.append(block)
                block = [layer_j]
        else:
            blocks.append(block)

    # Initialize
    T_12_block = []
    R_12_block = []
    T_21_block = []
    R_21_block = []
    
    n_ws = len(ws)    
    T_12_total = np.full((2, n_ws), 1, dtype=float)
    T_21_total = np.full((2, n_ws), 1, dtype=float)
    R_12_total = np.full((2, n_ws), 0, dtype=float)
    R_21_total = np.full((2, n_ws), 0, dtype=float)

    for block in blocks:        
        tr = get_tr_matrix(block, n_0, q_0, ws, sp=sp)
        T_12_block.append(tr[0])
        R_12_block.append(tr[1])
        T_21_block.append(tr[2])    
        R_21_block.append(tr[3])
    T_12_block = np.array(T_12_block)
    R_12_block = np.array(R_12_block)
    T_21_block = np.array(T_21_block)
    R_21_block = np.array(R_21_block)

    for i in range(len(blocks)):
        norm = 1 - R_21_total[0] * R_12_block[i]
        T_12_total[1] = T_12_total[0] * T_12_block[i] / norm
        R_12_total[1] = R_12_total[0] + T_12_total[0] * T_21_total[0] * R_12_block[i] / norm
        T_21_total[1] = T_21_total[0] * T_21_block[i] / norm
        R_21_total[1] = R_21_block[i] + T_12_block[i] * T_21_block[i] * R_21_total[0] / norm

        T_12_total[0] = T_12_total[1]
        R_12_total[0] = R_12_total[1]
        T_21_total[0] = T_21_total[1]
        R_21_total[0] = R_21_total[1]

    T_total = T_12_total[1]
    R_total = R_12_total[1]
    return {'T': T_total, 'R': R_total}