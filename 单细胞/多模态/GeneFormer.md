## 虚拟扰动

#### 1.背景

```shell
## 1.nlp 
## 2.bert 
## 深度学习模型用的都是较为基础的，不代码介绍
```

#### 2.基础（给单细胞分类）

```python 
##### 数据集就是单细胞分类的，
import torch

print(torch.cuda.is_available())

import torch
import numpy as np
from datasets import load_from_disk

dataset = load_from_disk('../data/hPBMC_disco.dataset/')

labels_train = dataset['train'].unique('cell_subtype')
labels_test = dataset['test'].unique('cell_subtype')
labels = list(set(labels_train) | set(labels_test))

label_dict = {labels[v]: v for v in range(len(labels))}

dataset = dataset.map(lambda x: {'label': label_dict[x['cell_subtype']]}, num_proc=4)


def collate_fn(data):
    X = [(np.array(i['input']) > 0).astype(float) for i in data]
    y = [i['label'] for i in data]
    X = torch.Tensor(np.array(X))  ## float32
    y = torch.LongTensor(y)  ## int64
    return X, y


data = dataset['train'].shuffle().select(range(10))

X, y = collate_fn(data)


class Model(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.fc = torch.nn.Sequential(
            torch.nn.Linear(in_features=4000, out_features=200), torch.nn.ReLU(), torch.nn.Dropout(p=0.4),
            torch.nn.Linear(in_features=200, out_features=100), torch.nn.ReLU(), torch.nn.Dropout(p=0.4),
            torch.nn.Linear(in_features=100, out_features=len(labels))
        )
    def forward(self, x):
        return self.fc(x)


model = Model() 

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

out = model(torch.randn(1, 4000))
out.argmax(dim=1)



loader = torch.utils.data.DataLoader(
    dataset['train'],
    batch_size=256,
    drop_last=True,
    shuffle=True,
    num_workers=10,
    collate_fn=collate_fn
)


def train(lr=1e-2, num_epochs=10, device=None):
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('training on', device)
    optimizer = torch.optim.SGD(model.parameters(), lr=lr)  
    criterion = torch.nn.CrossEntropyLoss()  
    model.train()
    if device != 'cpu':
        model.to(device)
    for epoch in range(num_epochs):
        for i, (X, y) in enumerate(loader):
            if device != 'cpu':
                X, y = X.to(device), y.to(device)
            out = model(X)
            loss = criterion(out, y)
            loss.backward()
            optimizer.step()

            optimizer.zero_grad()
        if epoch % 1 == 0:

            acc = (out.argmax(dim=1) == y).sum().item() / len(y)
            print(epoch, loss.item(), acc)
    return model

if __name__ == '__main__':

    model = Model()
    model = train(lr=0.1, num_epochs=5, device='cuda:0')

    torch.save(model.state_dict(), '../model/hPBMC_disco_celltype_classif.params')

    X, y = collate_fn(dataset['test'])

    model = Model()
    model.load_state_dict(torch.load('../model/hPBMC_disco_celltype_classif.params'))
    model.eval()

    with torch.no_grad():
        out = model(X)
        out = out.argmax(dim=1)
    correct = (out == y).sum().item()
    print(len(correct))
    print(y.shape)
    total = len(y[1])
    print(correct / total)


```
