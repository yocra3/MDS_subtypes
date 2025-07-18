import lightning as L
import torch
import sys
sys.path.append('./scripts/3.GeneticScores/pytorch_files/')

from GNNPatientSAGE import PatientGNNSAGE
from GNNPatientSAGE_v2 import PatientGNNSAGE_v2
from GNNPatientSAGE_v3 import PatientGNNSAGE_v3
from GNNPatientSAGE_v4 import PatientGNNSAGE_v4

class BaseEvalGNN(L.LightningModule):
    def __init__(self, ckpt_path, model_class, **model_kwargs):
        super().__init__()
        # init the pretrained LightningModule
        self.model = model_class(**model_kwargs)
        checkpoint = torch.load(ckpt_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])
        self.eval()
    
    def forward(self, x):
        return self.model(x)
    
    def predict(self, data):
        """Método unificado para hacer predicciones en grafos"""
        with torch.no_grad():
            return self(data).detach().numpy()

class EvalGNN(BaseEvalGNN):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__(ck_path, PatientGNNSAGE,
                        patient_feat_dim=patient_feat_dim, gene_feat_dim=gene_feat_dim,
                        hidden_gene_dim=hidden_gene_dim, hidden_dim=hidden_dim,
                        out_dim=out_dim, use_vaf=use_vaf, hidden_vaf=hidden_vaf,
                        learning_rate=learning_rate)

class EvalGNN2(BaseEvalGNN):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__(ck_path, PatientGNNSAGE_v2,
                        patient_feat_dim=patient_feat_dim, gene_feat_dim=gene_feat_dim,
                        hidden_gene_dim=hidden_gene_dim, hidden_dim=hidden_dim,
                        out_dim=out_dim, use_vaf=use_vaf, hidden_vaf=hidden_vaf,
                        learning_rate=learning_rate)

class EvalGNN3(BaseEvalGNN):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__(ck_path, PatientGNNSAGE_v3,
                        patient_feat_dim=patient_feat_dim, gene_feat_dim=gene_feat_dim,
                        hidden_gene_dim=hidden_gene_dim, hidden_dim=hidden_dim,
                        out_dim=out_dim, use_vaf=use_vaf, hidden_vaf=hidden_vaf,
                        learning_rate=learning_rate)

class EvalGNN4(BaseEvalGNN):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, 
                 dropout=0.3, l2_reg=1e-4, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__(ck_path, PatientGNNSAGE_v4,
                        patient_feat_dim=patient_feat_dim, gene_feat_dim=gene_feat_dim,
                        hidden_gene_dim=hidden_gene_dim, hidden_dim=hidden_dim,
                        out_dim=out_dim, dropout=dropout, l2_reg=l2_reg,
                        use_vaf=use_vaf, hidden_vaf=hidden_vaf,
                        learning_rate=learning_rate)