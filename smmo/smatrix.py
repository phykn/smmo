import numpy as np
from typing import List, Dict, Tuple, Any


class SMMO:
    '''Scattering Matrix method for Multilayer Optics'''
    
    def __init__(
        self,
        layers: List[Dict[str, Any]],
        config: Dict[str, Any]
    ) -> None:
        '''
        layers:
          - layer:
              n: refractive index
              k: absorption coefficient
              thickness: thickness of layer
              coherence: coherence in layer
        config:
          w: wavenumbers
          q: angle of incidence
          p: polarization
        '''
        if self.check_data(layers, config):
            if self.check_layers(layers):
                self.layers = layers
                self.n_0 = layers[0]["n"] + 1j * layers[0]["k"]
            if self.check_config(config):
                self.ws = config["w"]
                self.q_0 = config["q"]
                self.pol = config["p"]

    @staticmethod
    def check_data(
        layers: List[Dict[str, Any]],
        config: Dict[str, Any]
    ) -> bool:
        lengths = [len(config["w"])]
        for layer in layers:
            lengths.append(len(layer["n"]))
            lengths.append(len(layer["k"]))
        assert len(np.unique(lengths)) == 1, "size mismatch"
        return True

    @staticmethod
    def check_layers(
        layers: List[Dict[str, Any]]
    ) -> bool:
        assert len(layers) > 1, "The minimum number of layers is 2."
        return True

    @staticmethod
    def check_config(
        config: Dict[str, Any]
    ) -> bool:
        assert config["q"] >= 0 and config["q"] < 90, "q must be set between 0 and 90."
        assert config["p"] in ["s", "p"], 'polarization not in ["s", "p"]'
        return True

    @staticmethod
    def get_cos_qi(
        n_0: np.ndarray, 
        n_i: np.ndarray, 
        q_0: float
    ) -> np.ndarray:
        '''
        Note:
            Return cosine(theta) in ith layer
        Args:
            n_0: refractive index of 0th layer (complex)
            n_i: refractive index of ith layer (complex)
            q_0: angle of incidence
        Return:
            cosine values of ith layer (complex)
        '''      
        q_0 = np.full(
            shape=len(n_0), 
            fill_value=complex(q_0*np.pi/180, 0), 
            dtype=complex
        )
        return np.sqrt(1-np.square((n_0*np.sin(q_0)/n_i)))

    @staticmethod
    def get_kz(
        wavenumbers: np.ndarray,
        n_i: np.ndarray,
        cos_qi: np.ndarray        
    ) -> np.ndarray:
        '''
        Args:
            wavenumbers: unit is cm-1 (float)
            n_i: refractive index of ith layer (complex)
            cos_qi: cosine values of ith layer (complex)
            
        Return:
            wavevector in z direction (complex)
        '''
        return 2*np.pi*n_i*wavenumbers*cos_qi

    @staticmethod
    def get_fresnel_coeff_ij(
        n_i: np.ndarray, 
        n_j: np.ndarray, 
        cos_qi: np.ndarray, 
        cos_qj: np.ndarray,
        rt='r',
        sp='s'
    ) -> np.ndarray:
        '''
        Note:
            Calculate Fresnel coefficient between n_i and n_j
        Args:
            n_i: refractive index of ith layer (complex)
            n_j: refractive index of jth layer (complex)
            cos_qi: cosine values of ith layer (complex)
            cos_qj: cosine values of jth layer (complex)
            rt: "r" (refelction) or "t" (transmission)
            sp: "s" (s-polarization) or "p" (p-polarization) 
        Return:
            Fresnel coefficient (complex array): rs, rp, ts, tp
        '''        
        if (rt == "r") and (sp == "s"):   
            return (n_i*cos_qi-n_j*cos_qj)/(n_i*cos_qi+n_j*cos_qj)
        elif (rt == "r") and (sp == "p"):
            return (n_j*cos_qi-n_i*cos_qj)/(n_j*cos_qi+n_i*cos_qj)
        elif (rt == "t") and (sp == "s"):
            return 2*n_i*cos_qi/(n_i*cos_qi+n_j*cos_qj)
        elif (rt == "t") and (sp == "p"):
            return 2*n_i*cos_qi/(n_j*cos_qi+n_i*cos_qj)    
    
    def get_smatrix_components(
        self,
        layers: List[Dict[str, Any]]
    ) -> Tuple[np.ndarray]:
        '''Return scattering matrix components (S_11, S_12, S_21, S_22)'''
        # add a pseudo final layer
        layer_f = layers[-1].copy()
        layer_f["thickness"] = 0.0
        layer_f["coherence"] = True
        layers.append(layer_f)

        # get numbers
        num_layer = len(layers)
        num_data = len(self.ws)

        # initialize
        n = np.zeros((num_layer, num_data), dtype=complex)
        d = np.zeros((num_layer), dtype=complex)   
        cos_q = np.zeros((num_layer, num_data), dtype=complex)
        phase = np.zeros((num_layer, num_data), dtype=complex)    
        for i, layer in enumerate(layers):
            n[i] = layer["n"] + 1j * layer["k"]
            d[i] = layer["thickness"]
            cos_q[i] = self.get_cos_qi(self.n_0, n[i], self.q_0)
            phase[i] = np.exp(1j * self.get_kz(n[i], cos_q[i], self.ws) * d[i])
            
        # get scattering matrix components
        S_11 = np.full((num_data), 1, dtype=complex)
        S_12 = np.full((num_data), 0, dtype=complex)
        S_21 = np.full((num_data), 0, dtype=complex)
        S_22 = np.full((num_data), 1, dtype=complex)
        for i in range(num_layer-1):
            j = i + 1        
            t = self.get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], rt='t', sp=self.pol)
            r = self.get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], rt='r', sp=self.pol)

            I_11 = 1 / t
            I_12 = r / t
            I_21 = I_12
            I_22 = I_11

            S_11 = phase[i] * S_11 / (I_11 - phase[i] * S_12 * I_21)
            S_12 = ((phase[i] * S_12 * I_22 - I_12) * phase[j]) / (I_11 - phase[i] * S_12 * I_21)
            S_21 = S_22 * I_21 * S_11 + S_21
            S_22 = S_22 * I_21 * S_12 + S_22 * I_22 * phase[j]

        return S_11, S_12, S_21, S_22

    def get_tr_matrix_components(
        self,
        layers: List[Dict[str, Any]]
    ) -> Tuple[np.ndarray]:
        '''Return the inverse of a diagonal TR matrix components (T_12, R_12, T_21, R_21)'''        
        # forward
        S_11, S_12, S_21, S_22 = self.get_smatrix_components(layers)
        T_12 = (S_11 * np.conjugate(S_11)).real
        R_12 = (S_21 * np.conjugate(S_21)).real
        # backward
        S_11, S_12, S_21, S_22 = self.get_smatrix_components(layers[::-1])
        T_21 = (S_11 * np.conjugate(S_11)).real
        R_21 = (S_21 * np.conjugate(S_21)).real

        return T_12, R_12, T_21, R_21

    def __call__(
        self
    ) -> Dict[str, np.ndarray]:
        '''Return transmission and reflection'''

        pseudo_layers = []
        for i, layer in enumerate(self.layers):
            if layer["coherence"]:
                pseudo_layers.append(layer)
            else:
                pseudo_layer = dict(
                    n=layer["n"],
                    k=layer["k"],
                    thickness=0.0,
                    coherence=True
                )
                pseudo_layers.append(pseudo_layer)
                pseudo_layers.append(layer)
                pseudo_layers.append(pseudo_layer)

        blocks = []
        block = pseudo_layers[0:1]
        for i in range(len(pseudo_layers)):
            if i != len(pseudo_layers)-1:
                layer_i = pseudo_layers[i]
                layer_j = pseudo_layers[i+1]  
                if layer_i["coherence"] == layer_j["coherence"]:
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
        
        num_data = len(self.ws)    
        T_12_total = np.full((2, num_data), 1, dtype=float)
        T_21_total = np.full((2, num_data), 1, dtype=float)
        R_12_total = np.full((2, num_data), 0, dtype=float)
        R_21_total = np.full((2, num_data), 0, dtype=float)

        for block in blocks:        
            tr = self.get_tr_matrix_components(block)
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