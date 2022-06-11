import torch

trans_t = lambda t: torch.tensor([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,t],
    [0,0,0,1]
], dtype=torch.float)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
k = trans_t(5.0).to(device)
print(k)
