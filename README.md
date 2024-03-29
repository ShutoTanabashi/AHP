# 階層分析法プログラム

製作者:棚橋秀斗

## 注意

このファイル`README.md`自体は2024/2に作成しました．
ソフトウェアをコーディングした時期とは大きく異なります．

プログラムの作成自体は2019年初頭です．
当時の技術的な未熟さ故，拙いコーディングとなっていますがご容赦ください．

本プログラムを実行あるいは流用したことにより生じる一切の事象について製作者は責任を負いません。
自己責任でお願いいたします。

## 階層分析法について

> 階層分析法（かいそうぶんせきほう）は、意思決定における問題の分析において、人間の主観的判断とシステムアプローチとの両面からこれを決定する問題解決型の意思決定手法。
> AHP (Analytic Hierarchy Process) とも呼ばれる。
> ピッツバーグ大学のThomas L. Saatyが提唱した。  
> [Wikipedia 階層分析法](https://ja.wikipedia.org/wiki/階層分析法) より引用

多数の観点で評価される複数の選択肢の中から選択を行う場合に用いる手法です．
複数の選択肢のうちから選ばれた2個の選択肢の比較を順番に行うことで，それぞれの選択肢のスコアを計算します．

## 含まれるファイルについて

* `GUI_AHP.py` - AHP実行時のGUI表示に関するプログラム．
* `AHP.py` - AH Pの計算処理に関するプログラム

## 実行方法

プログラム実行時は

```sh
python GUI_AHP.py
```

のように`GUI_AHP.py`の方を実行してください．

1. GUIが立ち上がったら選択肢を順に入力してください。  
  （選択肢を日本語で入力すると最終出力画面で表示が豆腐になる可能性があるのでご注意ください）
2. 全ての選択肢を入力し終えたら下のボタンをクリックしてください。
3. 次に判断基準を順に入力してください。
4. 全ての選択肢を入力し終えたら下のボタンをクリックしてください。
5. 判断基準が2個表示されますので、どちらの判断基準がどの程度重要であるかボタンを押してください。
  判断基準の組合せ全てに対して実施します。  
  （この時に判断基準の重要度に大きな矛盾があるとAHPを実行できません。再度の入力を求められる場合があります）
6. 判断規準と2個の選択肢が表示されますので、その判断基準においてどちらがどの程度優れているかボタンを押してください。
  判断基準と選択肢の組合せ全てに対して実施します。  
7. AHPの計算結果が表示されます。
  グラフは各選択肢のスコアを表示しています。
  値が高い方が選択肢として好ましいと考えられます。
