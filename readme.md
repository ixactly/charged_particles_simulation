# ECR is GOD
git を clone しよう
- gitをinstall
    - git clone https://github.com/ixactly/charged_particles_simulation.git を実行
- めんどくさかったら、リンクから緑色の"clone or download"からフォルダをダウンロード

## 環境構築
- python3.6.5以上
- pip install ~で各種インストール
    - anacondaとかは使わない方がいい
- 必要なパッケージ群
    - matplotlib
    - numpy
    - scipy
        - 実行した時に〇〇がありません！って言われたらそれをinstallしてください
- pip install ~　でinstallできる ex.)pip install matplotlib
    - もしかしたらpip3(or pip3.6) install ~
        - python installの環境による
- 僕が使ってたwin機は環境構築されてるけどこのコードはwin機で動きませ〜ん
    - 相対パスの部分のコード変えればできるけど．

## 使い方
- macOSかlinuxで動かしてください
- 電位のバイナリを相対パスで保存してるのでCHARGED-PARTICLES-SIMULATIONのディレクトリから実行してください
- 各種パラメータを用意する
    - 引き出し電極とアインツェルレンズ
- パラメータ
    - Arになってます
    - 初期電流とかもデータと照らし合わせながら適当に弄ってみてください

## 引出電極
Electrode/electrode_rad0.pyまたElectrode/electrode_radx.pyを実行
標準出力に従って各種パラメータを入力してください

- 電極間距離
- 印加電圧
- (角度) 

制約上25mmまでしか離せません．おしまい．

## アインツェルレンズ
Einzel_Lens/EinzelLens_Potential.pyを実行
お好きなレンズ電圧を入力してください

## 電荷軌道
使いたい引出電極とアインツェルレンズのバイナリデータがElectrode/Electrode_dataとEinzel_Lens/Einzel_Lens_dataにあることを確認したら次に進んでください

### 全体像を見たい時
particles_trajectory.pyを実行してください
### アインツェルレンズ入射前を見たい時
trajectory_before_lens.pyを実行してください
グラフのタイトルにレンズ入射前のビーム径が表示されるはずです

## まとめ
疲れた．



