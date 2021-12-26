# CEvNS-UCAS Simulation Framwork

## 介绍
UCAS CEvNS 项目 Simulation 软件框架

## 软件架构

Base on ROOT v6.22/08 and Geant4 v10.3.3

## 使用说明
安装
```bash
git clone https://gitee.com/marry_tim/cevns-ucas.git
cd cevns-ucas
mkdir build
cmake ..
make
```


### Simulation:
1.  如何开始模拟：  
    ```bash
    cd build/source/0-Simulation
    ./starXP -[v] [Simucard.card]
    ```
    `-[v]`选择将开启可视化界面，`-[Simucard.card]`将使用用户指定的.card文件作为模拟的输入。
2.  Simucard文件设置(示例请见Simucard.card):  
    `Simucard.card`文件指定软件将以何种形式模拟何种信号

    #### GDML文件设置：  
        DetectorGDML: DetectorConstruction.xml

    #### 输出文件设置：
        OutputFile: xxx.root
        OutputLevel: 1

    #### 模拟方式设置：   
        1. 使用内置能谱模拟：
        GenEvents: 产生指定事例数  
        
        2. 使用SimpleGun模拟：
        GunType:Simple  
        PGName: e+/e-/gamma/p+/n
        PGEnergy:(MeV)
        PGPos:x,y,z
        PGPol:s1,s2,s3
        PGMomDir:px,py,pz
        
        3. 使用mac文件模拟：
        GPSMacFile:xxx.mac

        3. 使用CDEX 的ExGPS mac文件模拟：
        ExGPSMacFile:xxx.mac

    #### 模拟信号选择：  
        PGType:信号名称
        PGParameters:信号参数，以`,`隔开  
        
        PGN:中子本底抽样



