import torch

def tanh(x):
    return (torch.exp(x) - torch.exp(-x)) / (torch.exp(x) + torch.exp(-x))


x = torch.tensor(1.0, requires_grad=True)
y = tanh(x)
y.backward()
print(x.grad)


x_grad = x.grad

z = x_grad.to('cpu').detach().numpy().copy()

print(z)