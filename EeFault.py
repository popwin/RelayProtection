import numpy as np
import matplotlib.pyplot as plt

class EeFault(): 
    """
    由于发生故障时，采样得到的量并部满足正弦稳态情况
    因此 采用全波傅氏算法计算故障时基波的 实部 和 虚部
    值得注意的是，傅立叶算法同样是在计算 稳态正弦量
    也就是说，所求得的每一个复数是该频率下的 实部和虚部
    因此，当存在直流衰减分量时，该算法有极大误差，需要配合FIR滤波器（一阶差分）
    """
    def __init__(self,endT,f,fs,N):
        self.endT = endT
        self.f = f
        self.fs = fs
        self.N = N
        self.knum = 5
        self.omega = 2*np.pi*f
        self.alpha = 0
        self.Ialpha = self.alpha - 60*np.pi/180
        self.A0 = 3
        self.tor = 50 #整个电网在故障点断开后的 tor=0.02
        self.Uf = 35 #故障相电压 kV 有效值
        self.U2f = 10
        self.U3f = 5
        self.If = 10 #故障相电流 kA 有效值
        self.I2f = 5
        self.z1 = np.complex(0.59,4.52) #正序阻抗 全线路
        self.z0 = np.complex(1.77,13.5) #零序阻抗 全线路
        self.z1_set = np.abs(self.z1)*0.8 #一段阻抗保护定值
        return
    
    def testPlot(self,x,y,type=1):
        #init plot
        self.fig = plt.figure()
        self.ax =self.fig.add_subplot(111)
        self.ax.set_autoscaley_on(True)
        self.ax.grid()
        if type == 1:
            #line
            self.ax.plot(x,y,color='b')
            self.ax.relim()
            self.ax.autoscale_view()
            self.ax.set_ylim(-600,800)
            self.fig.canvas.draw()
        elif type == 2:
            self.ax.bar(x,y,1.8,color='b')
            self.fig.canvas.draw()
            
    def FaultSam(self):
        """
        模拟故障时采样得到的值
        """
        uFault_sam = [] #电压采样值
        iFault_sam = [] #电流采样值
        tdata = [] #采样时间点
        for item,t in enumerate(np.arange(0,self.endT,1./self.fs)):
            tdata.append(t)
            #叠加了不同频率的电压采样值
            sampTmpU = self.A0*np.exp(-self.tor*t) + \
                      self.Uf*1.414*np.sin(self.omega*t+self.alpha) + \
                      self.U2f*1.414*np.sin(2*self.omega*t+self.alpha+np.pi/3) + \
                      self.U3f*1.414*np.sin(3*self.omega*t+self.alpha+np.pi/2)
            #叠加了不同频率的电流采样值
            sampTmpI = self.A0*np.exp(-self.tor*t) + \
                      self.If*1.414*np.sin(self.omega*t+self.Ialpha) + \
                      self.I2f*1.414*np.sin(2*self.omega*t+self.Ialpha-np.pi/3) 
            #每一个时刻t采样保存一个值         
            uFault_sam.append(sampTmpU) 
            iFault_sam.append(sampTmpI)
            
        #self.testPlot(tdata,uFault_sam)
        return (tdata,uFault_sam,iFault_sam)
    
    def DealData(self,Npoints,type=1):
        
        kthFreq = np.arange(self.knum) # 需要计算的 频率个数
        Re = [0 for i in range(self.knum)] #用于存放不同频率下的 实部
        Img = [0 for i in range(self.knum)] #用于存放不同频率下的 实虚部
        Nsum = np.sum(Npoints[:self.N])
        N_evensum = np.sum([Npoints[even_index] for even_index in np.arange(self.N) if (even_index)%2==0])
        N_oddsum = Nsum-N_evensum
        if type == 1:
            #采用外循环为每一个点 内循环为每一中频率
            #为要模拟 采样的过程 即采样是一个点一点点得到的
            for index in np.arange(self.N):
                for k in kthFreq:
                    #采样得到一个点后 对所有谐波分量进行计算
                    Re[k] += (Npoints[index]*np.sin(k*index*2*np.pi/self.N))*2/self.N #用滤波后的值进行计算基频分量
                    Img[k] += (Npoints[index]*np.cos(k*index*2*np.pi/self.N))*2/self.N #如果想要观察滤波前，可将Yfilter替换为samPoint

            #将每一个谐波分量组合成复数
        elif type == 2:
            #改进全波傅立叶方法，利用一个周期既全波傅立叶计算，又滤除直流衰减分量
            for k in kthFreq:
                for index in np.arange(self.N):
                    
                    Re[k] += (Npoints[index]*np.sin(k*index*2*np.pi/self.N))*2/self.N 
                    Img[k] += (Npoints[index]*np.cos(k*index*2*np.pi/self.N))*2/self.N
                
                B = N_oddsum/N_evensum
                A = (1-B)*Nsum
                fenziRe = 2*A*B*np.sin(k*2*np.pi/self.N)
                fenziImg = 2*A*(1-B*np.cos(k*2*np.pi/self.N))
                fenmu = (1-2*B*np.cos(k*2*np.pi/self.N)+np.power(B,2))*self.N
                deltaRe = fenziRe/fenmu
                deltaImg = fenziImg/fenmu
                Re[k] = Re[k] - deltaRe
                Img[k] = Img[k] - deltaImg
            
            """
            fenziRe = 2*A*B*np.sin(kindex*2*np.pi/N)
                fenziImg = 2*A*(1-B*np.cos(kindex*2*np.pi/N))
                fenmu = (1-2*B*np.cos(kindex*2*np.pi/N)+np.power(B,2))*N
                deltaRe = fenziRe/fenmu
                deltaImg = fenziImg/fenmu
                Re[kindex] = Re[kindex] - deltaRe
                Img[kindex] = Img[kindex] -deltaImg
            
            """
                
            
        
        return (Re,Img)
    
    def FullFourier(self,samPoint):
        """
        全波傅氏算法的核心
                   N-1
              2    ----
       Urk = --- * \
              N    /    u(i)*sin(ki*2pi/N)   
                   ----
                   i=0
                   
                   N-1
              2    ----
       Uik = --- * \
              N    /    u(i)*sin(ki*2pi/N)
                   ----
                   i=0
        k：k次谐波
        i：第i个时刻
        u(i)：第i个采样点
        N：一个周期的采样个数
        Urk：第k次谐波的 实部
        Uik：第k次谐波的 虚部
        """        
          
        knFreqReIm = [] #其内每一个值对应着 该频率下 正弦量的 复数形式 
        
        
        #调用滤波函数
        #Fir = FIR() 
        #给定2N个采样值 前N个用于滤波 得到Yfilter大约是N个点
        #Yfilter = Fir.Casca3FIR(samPoint[:2*N])
        
        (Re,Img) = self.DealData(samPoint,type=2)
        
        for re,imag in zip(Re,Img):
            knFreqReIm.append(np.complex(re,imag))
        
        return knFreqReIm
    
    def FaultSamReIm(self,faultSample=None):
        """
        计算故障时 采样值中 不同频率下 正弦量的 实部和虚部
        有了故障时刻得到的 基频分量的 复数形式 
        接下来就可以进行 各种保护的 动作判据计算
        """
        
        #获得采样值
        if faultSample == None:
            (tdata,uFault_sam,iFault_sam) = self.FaultSam() #
        else:
            (uFault_sam,iFault_sam) = faultSample
        #求采样值中 不同频率 电压正弦量的 复数形式 该复数的模值为正弦量的 最大值
        UfnReIm = self.FullFourier(uFault_sam)
        #求采样值中 不同频率 电流正弦量的 复数形式
        IfnReIm = self.FullFourier(iFault_sam)
        
        #x = (np.arange(self.knum))*self.f
        #y = np.abs(IfnReIm)
        #self.testPlot(x,y,type=2)
        
        return (UfnReIm,IfnReIm)
    
    def FaultImped(self,u,i):
        """
        通过给定故障时计算得到的 电压和电流的复数，求故障时的阻抗
        该阻抗并非距离保护计算得到的阻抗，这里仅仅为简化而已
        Uma = Ua
        Ima = Ia + K*3I0
            Ur*Ir+Ui*Ii
        R = -----------
             Ir^2+Ii^2
             
            Ui*Ir-Ur*Ii
        Q = -----------
             Ir^2+Ii^2
        """
        #阻抗
        #Resis = (u.real*i.real+u.imag*i.imag)/(np.power(i.real,2)+np.power(i.imag,2))
        #React = (u.imag*i.real-u.real*i.imag)/(np.power(i.real,2)+np.power(i.imag,2))
        K = (self.z0-self.z1)/3/self.z1
        Um = u
        Im = i + K*i
        
        imped = Um/Im
        return imped
    
    def DistancePro(self,baseU,baseI):
        """
        距离保护
        """
        DistancePro_mark = False
        #计算故障时的基频分量阻抗值
        FaultImped = self.FaultImped(baseU,baseI)
        
        if np.abs(FaultImped) < self.z1_set:
            DistancePro_mark = True
        
        return (DistancePro_mark,FaultImped)
    
    
    
    
    
    