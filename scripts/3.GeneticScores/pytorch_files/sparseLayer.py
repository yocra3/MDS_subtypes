import torch
from torch import nn
import math

class SparseGO(torch.nn.Module):

    def __init__(self, input_genes, middle_layer, latent_go, relation_dict, bias = True, device = None, dtype = None):
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__()
        self.input_genes = input_genes
        self.middle_layer = middle_layer
        self.latent_go = latent_go
        self.relation_dict = relation_dict

        ## define multiplier
        nmult = middle_layer//latent_go

        ## create sparse weight matrix according to GO relationships
        mask = torch.zeros((self.input_genes, self.middle_layer), **factory_kwargs)
        for i in range(self.input_genes):
            for latent_go in self.relation_dict[i]:
                
                ## if we have middle layers larger than latent_go, we impose
                ## a layer size of a multiple of the latent space, and we will
                ## just duplicate the latent space
                for j in range(nmult):
                    mask[i, nmult*latent_go+j] = 1

        self.mask = mask
        self.weight = nn.Parameter(torch.empty((self.middle_layer, self.input_genes), **factory_kwargs))

        if bias:
             self.bias = nn.Parameter(torch.empty(middle_layer, **factory_kwargs))
        else:
             self.register_parameter('bias', None)

        self.reset_parameters()

    def forward(self, x):
        return (x @ ((self.weight * self.mask.T).T + self.bias))
        #return (torch.sparse.mm(x, self.weight))
        
    
    def reset_parameters(self) -> None:
        # Setting a=sqrt(5) in kaiming_uniform is the same as initializing with
        # uniform(-1/sqrt(in_features), 1/sqrt(in_features)). For details, see
        # https://github.com/pytorch/pytorch/issues/57109
        nn.init.kaiming_uniform_(self.weight, a=math.sqrt(5))
        if self.bias is not None:
            fan_in, _ = nn.init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1 / math.sqrt(fan_in) if fan_in > 0 else 0
            nn.init.uniform_(self.bias, -bound, bound)