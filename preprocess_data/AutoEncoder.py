import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

from torch.utils.data import DataLoader

import matplotlib.pyplot as plt

import time
import os
import copy
import math
import glob

from collections import deque
from pathlib import Path


#-----load data-----#
EXP_df = pd.read_csv("C:/Users/tomsh/OneDrive/デスクトップ/研究テーマ/TCGAPAADdata/preprocess_data/PAAD_Exp.csv")


#-----preprocess function-----#
class AEDataset():
    def __init__(self, X, y):
        self.X = torch.Tensor(X.values)
        self.y = torch.Tensor(y.values)

    def __len__(self):
        return len(self.y)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx] 
        #後でdata_loader作るときに、シャッフルした方を"train"にして、シャッフルしてないほうを"test"にして辞書にする
        #キーからデータを取ってこれるようにgetitem関数作る


#-----Autoencoder-----#
class Encoder(nn.Module):
    def __init__(self, input_features, encoding_dim): #encoding_dimはリストで渡す
        super().__init__()
        self.fc1 = nn.Linear(input_features, encoding_dim[0])
        self.pool1 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(encoding_dim[0], encoding_dim[1])
        self.pool2 = nn.Dropout(0.5)
        print

    def forward(self, x):
        #print("Encoder入った")
        #print("Encoder入力のxのsizeは:{}".format(x.shape))
        x = torch.tanh(self.fc1(x)) #ReLUだったら入力負の時0になっちゃうからいったんtanhでやってみる。元の論文の原因これReLUにしてたからでは？
        x = self.pool1(x)
        x = torch.tanh(self.fc2(x))
        x = self.pool2(x)
        #print("Encoder出力のxのsizeは:{}".format(x.shape))
        #print("Encoder出た")
        return x
    
class Decoder(nn.Module):
    def __init__(self, encoding_dim, input_features):
        super().__init__()
        self.fc3 = nn.Linear(encoding_dim[1], encoding_dim[0])
        self.fc4 = nn.Linear(encoding_dim[0], input_features)

    def forward(self, x):
        x = torch.tanh(self.fc3(x))
        x = torch.tanh(self.fc4(x))
        return x
    
class AutoEncoder(nn.Module):
    def __init__(self, input_features, encoding_dim):
        super().__init__()
        self.encoder = Encoder(input_features, encoding_dim)
        self.decoder = Decoder(encoding_dim, input_features)

    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x
    

#-----train関数-----#
def train_model(model, loss_func, optimizer, data_loader, n_epochs, fout, device):
    #初期化
    pkl_queue = deque()
    best_loss = 100.0
    best_epoch = 0
    best_model_weights = model.state_dict()  # state_dict()はtorchの関数
    since = time.time()
    end = time.time()

    print(model, "\n")

    for epoch in range(n_epochs):
        print("EPOCH:{}/{}".format(epoch+1, n_epochs), end="")
        print("EPOCH:{}/{}".format(epoch+1, n_epochs), end="", file=fout)

        for phase in ["train"]:
            model.train(True)

            #データの指定
            data = data_loader[phase]

            #初期化
            running_loss = 0

            #ミニバッチに対するループ処理
            for _, (data_train, target_train) in enumerate(data):
                optimizer.zero_grad()
                x = data_train.to(device)
                y = target_train.to(device)

                with torch.set_grad_enabled(phase == "train"):
                    y_pred = model(x)
                    loss = loss_func(y_pred, y)

                    loss.backward()
                    optimizer.step()

                running_loss += loss.item()
                
            epoch_loss = running_loss / (len(data)/len(x))

            #最も損失が小さかったモデルを保存
            if epoch_loss < best_loss:
                best_loss = epoch_loss
                best_epoch = epoch
                best_model_weights = copy.deepcopy(model.state_dict())
                torch.save(best_model_weights, "{}_epoch{}.pkl".format(fout.name.split(".txt")[0], epoch+1))
                pkl_queue.append("{}_epoch{}.pkl".format(fout.name.split(".txt")[0], epoch+1))

                if len(pkl_queue) > 1:
                    pkl_file = pkl_queue.popleft()
                    os.remove(pkl_file)
            
            #予測の出力
            print(",{}Loss:{:.4f}, Time:{:.4f}".format(phase, epoch_loss, time.time()-end), end="")
            print(",{}Loss:{:.4f}, Time:{:.4f}".format(phase, epoch_loss, time.time()-end), end="", file=fout)
            print("\n", end="")
            print("\n", end="", file=fout)

            end = time.time()

    #トレーニング結果の表示
    time_elapsed = time.time() - since
    print("\nTraining completed in {:.0f}m {:.0f}s".format(time_elapsed // 60, time_elapsed % 60))
    print("Best loss: {:.4f} at epoch {}".format(best_loss, best_epoch))


#-----Epoch-lossグラフの描画-----#
def plot_loss(file, name):

    #トレーニングログファイルの読み込み
    df = pd.read_csv(file, header=None, sep=r"\s+")
    for col in df.columns:
        df[col] = df[col].str.replace('[A-Za-z]+:', '', regex=True)
        df[col] = df[col].str.replace(',', '', regex=True)    

    #線グラフを作成
    fig, ax = plt.subplots()
    plt.plot(range(len(df)), df.iloc[:, 1].astype('float'))
    ax.set_title(f"MSE loss for \n{file}")
    ax.set_xlabel("EPOCHS")
    ax.set_ylabel("MSE loss")
    fig.savefig(f"{name}.png")
    plt.close(fig)


#-----Eval model-----#
def eval_model(model, data_loader, device):
    #初期化
    running_mse = 0
    preds = []

    data = data_loader["test"]
    model.eval()

    for data_test, target_test in data:
        x = data_test.to(device)
        y = target_test.to(device)

        #予測
        with torch.no_grad():
            y_pred = model(x)

            #MSEの計算
            sq_loss = ((y_pred - y)*(y_pred - y)).sum().data
            running_mse += sq_loss

            preds.append(y_pred[0])
    #予測スコアを表示
    preds = np.vstack(preds)
    mse = math.sqrt(running_mse / len(data))
    print("MSE: {}".format(mse))

    return preds


#-----Hyper parameters-----#
batch_size = 4
n_epochs = 10
lr = 0.0008


#-----Check file path-----#
def check_path(filename):
    if not Path(filename).is_dir():
        Path(filename).parents[0].mkdir(parents=True, exist_ok=True)
    return filename


#-----Main function-----#
def main():
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    path = Path.cwd()

    #input_dataの読み込み
    df = pd.read_csv(Path(path, "PAAD_Exp.csv"))
    X = df.astype(np.float32)

    #load AutoEncoder
    model = AutoEncoder(3000, [200, 50])
    model = model.to(device)

    #dataset 作る
    AEdata = AEDataset(X, X)
    train_data = DataLoader(AEdata, batch_size=batch_size, shuffle=True)
    test_data = DataLoader(AEdata, batch_size=batch_size, shuffle=False)
    data_loader = {"train": train_data, "test" : test_data}

    loss_func = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    #トレーニングログファイルの設定
    train_log = "AE_lr{}_epochs{}_batch{}.txt".format(lr, n_epochs, batch_size)
    fout = open(check_path(str(Path(path, "AutoEncoder", train_log))), "w") #ログファイルを書き込みモードで開く

    #train AutoEncoder
    train_model(model, loss_func, optimizer, data_loader, n_epochs, fout, device)

    fout.close() #書き込み終了

    #train lossの書き出し
    plot_loss(str(Path(path, "AutoEncoder", train_log)), check_path(str(Path(path, "AutoEncoder", f"epoch{n_epochs}_loss"))))

    #最も性能の良いAEモデルの読み込み
    trained_model = glob.glob(str(Path(path, "AutoEncoder", "{}_*.pkl".format(train_log.split(".txt")[0]))))[0]
    model.load_state_dict(torch.load(trained_model), strict=False)

    #学習済みモデルを用いた予測
    decoded_result = eval_model(model, data_loader, device)
    print(decoded_result)

    #もっともよい学習済みモデルからの特徴抽出
    bottleneck_features = model.encoder(torch.Tensor(X.values)).detach().numpy()
    print(bottleneck_features)

    #低次元特徴量の保存
    np.savetxt(check_path(str(Path(path, "Bottleneck", f"epoch{n_epochs}_std_Exp.csv"))), bottleneck_features, delimiter=",")
    np.save(check_path(str(Path(path, "Bottleneck", f"epoch{n_epochs}_std_Exp.npy"))), bottleneck_features)

    #decodeしたデータの保存
    np.savetxt(check_path(str(Path(path, "DecodedExp", f"epoch{n_epochs}_std_Exp.csv"))), decoded_result, delimiter=",")


if __name__ == "__main__":
    main()
    
