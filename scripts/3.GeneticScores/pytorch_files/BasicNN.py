from torch import nn
from pycox.models.loss import cox_ph_loss
from lifelines.utils import concordance_index
from torch.optim import Adam
import lightning as L

# Define model
class BasicNN(L.LightningModule):
    def __init__(self, inp_dim, hidden_dim, learning_rate):
        super().__init__()
        self.layer_1 = nn.Sequential(nn.Linear(inp_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))
        self.predic = nn.Sequential(nn.Linear(hidden_dim, 1))
        self.learning_rate = learning_rate
    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.layer_1(x)
        x_hat = self.predic(z)
        loss = cox_ph_loss(x_hat, y[:, 0], y[:, 1])
        c_index = concordance_index(y[:, 0].cpu().numpy(), -x_hat.cpu().detach().numpy(), y[:, 1].cpu().numpy())
        self.log("train_loss", loss, prog_bar=True, on_epoch=True)
        self.log("train_c_index", c_index, prog_bar=True, on_epoch=True, on_step = False)
        return loss
    def validation_step(self, batch, batch_idx):
        # this is the validation loop
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.layer_1(x)
        x_hat = self.predic(z)
        val_loss = cox_ph_loss(x_hat, y[:, 0], y[:, 1])
        c_index = concordance_index(y[:, 0].cpu().numpy(), -x_hat.cpu().detach().numpy(), y[:, 1].cpu().numpy())
        self.log("val_loss", val_loss, prog_bar=True)
        self.log("val_c_index", c_index, prog_bar=True)
    def configure_optimizers(self):
        optimizer = Adam(self.parameters(), lr=self.learning_rate)
        return optimizer
    def forward(self, x):
        x = x.view(x.size(0), -1)
        z = self.layer_1(x)
        x_hat = self.predic(z)
        return x_hat


