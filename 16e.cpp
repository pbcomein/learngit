// 802.16E.cpp : 定义控制台应用程序的入口点。
//git第三次修改
// ProtocolEOOA.cpp : implementation file
#include "stdafx.h"
#include "ProtocolE.h"
#include "SignalSource.h"
#include "filter_para.h"
#include <string>
#include <iostream>
#include <fstream>
#include "fft\vectormath.h"
#include "fft\fft.h"

using namespace splab;
ProtocolE CProtocol;

parameters_16e::~parameters_16e()
 {
 }

 ProtocolE::ProtocolE()
{
	QPSK[0]=complex<double>(0.7071, 0.7071);
	QPSK[1] = complex<double>(0.7071, -0.7071);
	QPSK[2] = complex<double>(-0.7071, +0.7071);
	QPSK[3] = complex<double>(-0.7071, -0.7071);
	_16QAM[0]=complex<double>(0.3162, 0.3162); _16QAM[1]=complex<double>(0.3162, 0.9487);
	_16QAM[2]=complex<double>(0.3162, -0.3162);  _16QAM[3]=complex<double>(0.3162, -0.9487);
	_16QAM[4]=complex<double>(0.9487, 0.3162); _16QAM[5]=complex<double>(0.9487, 0.9487);
	_16QAM[6]=complex<double>(0.9487, -0.3162);  _16QAM[7]=complex<double>(0.9487, -0.9487);
	_16QAM[8]=complex<double>(-0.3162, 0.3162);  _16QAM[9]=complex<double>(-0.3162, 0.9487);
	_16QAM[10]=complex<double>(-0.3162, -0.3162);  _16QAM[11]=complex<double>(-0.3162, -0.9487);
	_16QAM[12]=complex<double>(-0.9487, 0.3162); _16QAM[13]=complex<double>(-0.9487, 0.9487);
	_16QAM[14]=complex<double>(-0.9487, -0.3162);  _16QAM[15]=complex<double>(-0.9487, -0.9487);

	_64QAM[0]=complex<double>(0.4629, 0.4629); _64QAM[1]=complex<double>(0.4629, 0.1543);
	_64QAM[2]=complex<double>(0.4629, 0.7715); _64QAM[3]=complex<double>(0.4629, 1.0801);
	_64QAM[4]=complex<double>(0.4629, -0.4629);  _64QAM[5]=complex<double>(0.4629, -0.1543);
	_64QAM[6]=complex<double>(0.4629, -0.7715);  _64QAM[7]=complex<double>(0.4629, -1.0801);

	_64QAM[8]=complex<double>(0.1543, 0.4629); _64QAM[9]=complex<double>(0.1543, 0.1543);
	_64QAM[10]=complex<double>(0.1543, 0.7715); _64QAM[11]=complex<double>(0.1543, 1.0801);
	_64QAM[12]=complex<double>(0.1543, -0.4629);  _64QAM[13]=complex<double>(0.1543, -0.1543);
	_64QAM[14]=complex<double>(0.1543, -0.7715);  _64QAM[15]=complex<double>(0.1543, -1.0801);

	_64QAM[16]=complex<double>(0.7715, 0.4629); _64QAM[17]=complex<double>(0.7715, 0.1543);
	_64QAM[18]=complex<double>(0.7715, 0.7715); _64QAM[19]=complex<double>(0.7715, 1.0801);
	_64QAM[20]=complex<double>(0.7715, -0.4629);  _64QAM[21]=complex<double>(0.7715, -0.1543);
	_64QAM[22]=complex<double>(0.7715, -0.7715);  _64QAM[23]=complex<double>(0.7715, -1.0801);

	_64QAM[24]=complex<double>(1.0801, 0.4629); _64QAM[25]=complex<double>(1.0801, 0.1543);
	_64QAM[26]=complex<double>(1.0801, 0.7715); _64QAM[27]=complex<double>(1.0801, 1.0801);
	_64QAM[28]=complex<double>(1.0801, -0.4629);  _64QAM[29]=complex<double>(1.0801, -0.1543);
	_64QAM[30]=complex<double>(1.0801, -0.7715);  _64QAM[31]=complex<double>(1.0801, -1.0801);

	_64QAM[32]=complex<double>(-0.4629, 0.4629); _64QAM[33]=complex<double>(-0.4629, 0.1543);
	_64QAM[34]=complex<double>(-0.4629, 0.7715); _64QAM[35]=complex<double>(-0.4629, 1.0801);
	_64QAM[36]=complex<double>(-0.4629, -0.4629);  _64QAM[37]=complex<double>(-0.4629, -0.1543);
	_64QAM[38]=complex<double>(-0.4629, -0.7715);  _64QAM[39]=complex<double>(-0.4629, -1.0801);

	_64QAM[40]=complex<double>(-0.1543, 0.4629); _64QAM[41]=complex<double>(-0.1543, 0.1543);
	_64QAM[42]=complex<double>(-0.1543, 0.7715); _64QAM[43]=complex<double>(-0.1543, 1.0801);
	_64QAM[44]=complex<double>(-0.1543, -0.4629);  _64QAM[45]=complex<double>(-0.1543, -0.1543);
	_64QAM[46]=complex<double>(-0.1543, -0.7715);  _64QAM[47]=complex<double>(-0.1543, -1.0801);

	_64QAM[48]=complex<double>(-0.7715, 0.4629); _64QAM[49]=complex<double>(-0.7715, 0.1543);
	_64QAM[50]=complex<double>(-0.7715, 0.7715); _64QAM[51]=complex<double>(-0.7715, 1.0801);
	_64QAM[52]=complex<double>(-0.7715, -0.4629);  _64QAM[53]=complex<double>(-0.7715, -0.1543);
	_64QAM[54]=complex<double>(-0.7715, -0.7715);  _64QAM[55]=complex<double>(-0.7715, -1.0801);

	_64QAM[56]=complex<double>(-1.0801, 0.4629); _64QAM[57]=complex<double>(-1.0801, 0.1543);
	_64QAM[58]=complex<double>(-1.0801, 0.7715); _64QAM[59]=complex<double>(-1.0801, 1.0801);
	_64QAM[60]=complex<double>(-1.0801, -0.4629);  _64QAM[61]=complex<double>(-1.0801, -0.1543);
	_64QAM[62]=complex<double>(-1.0801, -0.7715);  _64QAM[63]=complex<double>(-1.0801, -1.0801);


}
void ProtocolE::tx_param()
{
/*******************************************************************************************/
/**                                  System parameters                                    **/
/*******************************************************************************************/
//OFDM基本参数
    para.NumOfUsers = 1;
	
	para.n_factor = (float)28/25; //抽样因子

    para.Fs = floor(para.n_factor * para.BW / 8000) * 8000; //抽样频率
    para.Subcarrier_spacing = para.Fs / para.FFT_Size; //子载波频率间隔
    para.Tb = 1 / para.Subcarrier_spacing; //OFDM有效符号持续时间/子载波频率间隔为有效符号持续时间的倒数
    para.Tg = para.Tb * para.G_factor; //循环前缀的持续时间
    para.Ts = para.Tb + para.Tg; //OFDMA符号时间
    para.Sampling_Time = para.Tb / para.FFT_Size; //信号时域抽样间隔T/N
    para.LengthOfTgSp = para.G_factor * para.FFT_Size; //循环前缀的长度
    para.LengthOfTsSp = para.FFT_Size + para.LengthOfTgSp; //OFDMA符号的长度
    para.SymbolsForPreamble = 1; //前导码所占的OFDMA符号数
	para.SymbolsPerFrame = para.SymbolsForPreamble+para.SymbolsForBurst; //下行子帧符号总数

	para.NumPerSlot = para.LengthPerSubchannel * 2; //每个时隙子信道数量。对于下行PUSC 一个时隙占用一个子信道两个OFDM符号时间
	para.NumOfUnit = para.LengthPerSubchannel * para.Usedsubchannel; //数据单元个数=数据子载波个数（尝试添加子信道选择）
    
    para.DurationPerFrame = para.Ts * para.SymbolsPerFrame; //下行帧数据时长
	
    para.LengthOfFrame = para.LengthOfTsSp * para.SymbolsPerFrame; //每帧数据位长


/********************************************************************************************
                                信道编码参数的初始化
********************************************************************************************/
	switch(para.MCS_index){
	case 0:
		para.CodeRate = 0.5;para.Modulation_Mode=4; para.j_block=6;break;
	case 1:
		para.CodeRate = 0.75;para.Modulation_Mode=4;para.j_block=4;break;
	case 2:
		para.CodeRate = 0.5;para.Modulation_Mode=16;para.j_block=3;break;
	case 3:
		para.CodeRate = 0.75;para.Modulation_Mode=16;para.j_block=2;break;
	case 4:
		para.CodeRate = 0.5;para.Modulation_Mode=64;para.j_block=2;break;
	case 5:
		para.CodeRate = 2/3.0;para.Modulation_Mode=64;para.j_block=1;break;
	case 6:
		para.CodeRate = 0.75;para.Modulation_Mode=64;para.j_block=1;break;
	case 7:
		para.CodeRate = 0.5;para.Modulation_Mode=4;para.j_block=6;break;
	case 8:
		para.CodeRate = 0.75;para.Modulation_Mode=4;para.j_block=4;break;
	case 9:
		para.CodeRate = 0.5;para.Modulation_Mode=16;para.j_block=3;break;
	case 10:
		para.CodeRate = 0.75;para.Modulation_Mode=16;para.j_block=2;break;
	case 11:
		para.CodeRate = 0.5;para.Modulation_Mode=64;para.j_block=2;break;
	case 12:
		para.CodeRate = 2/3.0;para.Modulation_Mode=64;para.j_block=1;break;
	}

    para.CBSize = 48 * (float)log((float)para.Modulation_Mode) / (float)log((float)2); //NumOfSolt每个时隙占用的子信道数，NumPerSlot每个时隙的调制后数据个数（48），信道编码后数据块比特数
    para.UncBSize = (int)(para.CBSize * para.CodeRate); //信道编码前数据块比特数

	Preamble_generate(para.Preamble_out);//产生前导
	poilts_generate(para.Poilts_Out);
/********************************************************************************************/
/**                              FCH DL_MAP UL_MAP 生成相关参数                           **/
/**                                      DL_MAP 144                                        **/
/**                                      UL_MAP 48                                         **/
/**                                      FCH    192                                        **/
/********************************************************************************************/
    /*********************************************/
    /**             DL_MAP生成相关参数         **/
    /*********************************************/
    para.Management_Message_Type = 2;                     //默认
    para.Frame_Duration_Code = 4;    //帧时长5ms
    para.Frame_Number_index = 1;                          //默认
    para.BS_ID_index = 1;
    para.NumOfDL_Frame_Symbol = para.SymbolsPerFrame;//下行子帧OFDMA符号数

    /**********************************************
                        DL_MAP_IE
    ***********************************************/
    para.DIUC_index = para.MCS_index; //0-12
    para.OFDMA_Symbol_offset = 3; //第二个突发符号偏移，第一个突发默认为1
//    para.Subchannel_offset_burst2 = 0;//第二个突发子信道偏移

	int Subchannel_offset0=9;  //第一个突发子信道偏移
	if(para.DL_MAP_repetition==0)
	{
		Subchannel_offset0=4+5;
	}
		
	else
	{
		Subchannel_offset0=4+5*para.DL_MAP_repetition;
	}
		

    if (para.DL_MAP_Source_mode)
    {
		para.DL_MAP = FromFile(100,para.DL_MAP_Source_filename);
    }
	else
	{
		DL_MAP_generate(para.Power_boost, para.DIUC_index, para.OFDMA_Symbol_offset,Subchannel_offset0, para.Subchannel_offset_burst2, para.SymbolsForBurst2, para.Subchannel_num_burst2,para.Burst_repetition, para.DL_MAP);
	}

	MAC_head(0,38,2,para.MAC_head_DLMAP);

	vector<int>DL_MAP_PDU;
	DL_MAP_PDU.insert(DL_MAP_PDU.end(),para.MAC_head_DLMAP.begin(),para.MAC_head_DLMAP.end());
	DL_MAP_PDU.insert(DL_MAP_PDU.end(),para.DL_MAP.begin(),para.DL_MAP.end());
    DL_MAP_modulation(DL_MAP_PDU, para.DL_MAP_repetition, para.DL_MAP_Coding, para.DL_MAP_out);// DL_MAP_modulation(); //144

	
    /**********************************************
              FCH生成相关参数及生成结果
    ***********************************************/

    para.DL_MAP_Length = para.DL_MAP_out.size()/ 48; //以时隙为单位*/

    FCH_generate(para.DL_Frame_Prefix,para.Group);
    FCH_modulation(para.DL_Frame_Prefix, para.FCH_out); //192

/********************************************************************************************
                                Subcarrier Allocation
********************************************************************************************/
    PUSC();
	


/********************************************************************************************
                     MIMO Multipath_Channel Parameters and initialization
********************************************************************************************/
    para.NumofFreq = 8;
    para.NumofChannels = para.Nt * para.Nr; //MIMO信道个数
    para.Channel_Type;//选择信道模式
    para.uplink = 0; //Downlink
    para.MobileSpeed = 3; //km/h
    para.SampleIndex = 0;
    for (int i = 0; i < para.Nt * para.Nr; i++)
    {
        para.Next_SampleIndex.push_back(0);
    }
    para.corrmodel = 0;
    para.Data_From_Pre.push_back(complex <float> (0, 0.00000000001));
    para.UpdateInterval = para.LengthOfFrame;
    para.UpdatesPerBurst = 1;
	//信道初始化输出参数初始化
	para.UpdateInterval;
	para.UpdatesPerBurst;
	para.dt;
	para.NumOfTaps;
    para.aglsprR;
	para.aglsprT;
	para.dR;
	para.dT;
	para.corrmodel;
	para.alf;
	para.fd;
	para.loson;
}

void ProtocolE::Run( vector<complex<double>>  &Data_Tx,vector<complex<double>>  &Data_Tx_2 )//生成Data_Tx和Data_Tx_2,分别代表根据协议加cp后的正复数数据，然后传送至MainFrm画图，若增加信道，将信道输出后的信号输出。
{
	tx_param();
			/**********************************************
								Initialization
			***********************************************/
            vector <int> Data_In,burst1,burst2,burst3; //生成数据存储
            vector <int> InterLeavData1,InterLeavData2,InterLeavData3; //扰码、CC、交织后的生成数据存储
            vector <complex <double> > Mod_Out, Mod_Out_temp1,Mod_Out_temp2,Mod_Out_temp3; //
			vector <complex <double> > Rep_Out; //重复编码数据存储
            vector <complex <double> > Mod_Out_Encoded,Mod_Out_Encoded_2; //STBC antenna_1  STBC antenna_2
            vector <complex <double> > IFFT_DataIn,IFFT_DataIn_2; //OFDMA符号存储

		    /**********************************************
                            Adding the FCH DL_MAP
            DL_MAP的位置紧随FCH之后，用于向所有的用户终端广播下行子帧的资源分配情况
               
            ***********************************************/

            vector <complex <double> > Prefix;
			int num_of_prefix;
 
			Prefix.insert(Prefix.end(), para.FCH_out.begin(), para.FCH_out.end());
			if (para.DL_MAP_state)Prefix.insert(Prefix.end(), para.DL_MAP_out.begin(), para.DL_MAP_out.end());
			num_of_prefix = Prefix.size();		


            /**********************************************
                        Generate the input data
            ***********************************************/
			int real_downbits;
			vector<int> SignalSource;
			real_downbits=para.NumOfUnit * para.Subchannel_num_burst2/para.Usedsubchannel*(float)log((float)para.Modulation_Mode) / (float)log((float)2)*para.CodeRate;
				//为重复编码留空间

			//第一个突发数据
			SignalSource = GenerateSource(para.SignalSource_mode, para.NumOfUnit*2-num_of_prefix, para.SignalSource_seed, para.SignalSource_filename);
			for(int j=0;j<para.NumOfUnit*2-num_of_prefix;j++)
				burst1.push_back(SignalSource[j]);
			//第二个突发数据
            for (int i = 0; i < para.SymbolsForBurst2; i++)
            {
				SignalSource = GenerateSource(para.SignalSource_mode, real_downbits, para.SignalSource_seed, para.SignalSource_filename);
				for(int j=0;j<real_downbits;j++)
					burst2.push_back(SignalSource[j]);
				SignalSource.clear();
            }
			//第三个突发，即空时编码
			int real_downbits_burst3=para.NumOfUnit;

			 for (int i = 0; i <10; i++)
            {
				SignalSource = GenerateSource(para.SignalSource_mode, real_downbits_burst3, para.SignalSource_seed, para.SignalSource_filename);
				for(int j=0;j<real_downbits;j++)
					burst3.push_back(SignalSource[j]);
				SignalSource.clear();
            }

			block_channel_code(48,burst1,6,0,InterLeavData1);
			modulation(InterLeavData1,0,Mod_Out_temp1);

			
            block_channel_code(para.UncBSize, burst2,para.j_block,para.MCS_index,InterLeavData2);

			modulation(InterLeavData2, para.MCS_index,Mod_Out_temp2);

			repetition(para.Burst_repetition,Mod_Out_temp2,Rep_Out);

			block_channel_code(48,burst3,6,0,InterLeavData3);

			modulation(InterLeavData3, 0,Mod_Out_temp3);

			/**********************************************
                        UL_MAP&DL_MAP Inserting
            ***********************************************/
		
			Mod_Out.insert(Mod_Out.end(), Prefix.begin(), Prefix.end());
			Mod_Out.insert(Mod_Out.end(), Mod_Out_temp1.begin(),Mod_Out_temp1.end());
			Mod_Out.insert(Mod_Out.end(), Rep_Out.begin(),Rep_Out.end());

            STBC_encode(Mod_Out, para.Nt, para.NumOfUnit, Mod_Out_Encoded,Mod_Out_Encoded_2);//空时编码

			IFFT_DataIn.resize(para.FFT_Size * para.SymbolsForBurst); //初始化OFDMA符号结构

			Subcarrier_Mapping(Mod_Out_Encoded,IFFT_DataIn); //映射数据和导频



				vector <complex <double> > IFFT_DataIn_All_temp,IFFT_DataIn_All,IFFT_DataIn_All_2;

				IFFT_DataIn_All_temp.insert(IFFT_DataIn_All_temp.end(), para.Preamble_out.begin(), para.Preamble_out.end());
				IFFT_DataIn_All_temp.insert(IFFT_DataIn_All_temp.end(), IFFT_DataIn.begin(), IFFT_DataIn.end());
				

				int symbol_all = IFFT_DataIn_All_temp.size() / para.FFT_Size;

				for (int i = 0; i  < symbol_all; i++)
				{
					for (int j = para.FFT_Size/2; j < para.FFT_Size; j++) 
					{
						IFFT_DataIn_All.push_back(IFFT_DataIn_All_temp[i * para.FFT_Size + j]);
					}
					for (int j = 0; j < para.FFT_Size/2; j++)
					{
						IFFT_DataIn_All.push_back(IFFT_DataIn_All_temp[i * para.FFT_Size + j]);
					}
				}

				if(para.Nt== 2)
				{
					IFFT_DataIn_All_2.insert(IFFT_DataIn_All_2.end(), para.Preamble_out.begin(), para.Preamble_out.end());
				    IFFT_DataIn_All_2.insert(IFFT_DataIn_All_2.end(), IFFT_DataIn_2.begin(), IFFT_DataIn_2.end());
					if (para.SYNC_State == true) IFFT_DataIn_All_2.insert(IFFT_DataIn_All_2.end(), para.SYNC_out.begin(), para.SYNC_out.end());
					for (int i=IFFT_DataIn_All_2.size();i<=para.FFT_Size*14;i++)
					{
						IFFT_DataIn_All_2.push_back(0);
					}
				}
				/**********************************************
                                  CDD precoding
                ***********************************************/
				if(para.CDD_State)
				{
					 FreOffset_Add( IFFT_DataIn_All_2,para.CDD_Num,para.Sampling_Time);
				}

				IFFT_CP(IFFT_DataIn_All,Data_Tx);

				if(para.Nt== 2)
				{
					IFFT_CP(IFFT_DataIn_All_2,Data_Tx_2);
				}

				/**********************************************
                            Frequency Offset Adding
                ***********************************************/
				FreOffset_Add( Data_Tx,para.FrequencyOffset,para.Sampling_Time);
				if (para.Nt== 2)
					FreOffset_Add( Data_Tx_2,para.FrequencyOffset,para.Sampling_Time);
			    /**********************************************
                            Multipath Fading Channel
                ***********************************************/
				if (para.Channel_State)
				{
				vector <float> Path_Delay, Path_Average_Atmp, aglT, An, Bn, Wn, RndPhase, aglR;
				vector <int> distypeR, distypeT;
				//complex<double> ch_Data_Out(0,0); 

				////用于测试信道增益问题
				//para.LengthOfFrame = 1000;
				//int BW=20e+6;
				//para.DurationPerFrame= 5e-5;
				//para.SampleIndex = 0;
				//para.NumofFreq = 8;
				//para.uplink = 1;
				
				channle_init_996((float)para.LengthOfFrame, para.DurationPerFrame, para.MobileSpeed, para.NumofFreq, para.NumofChannels, para.Channel_Type, para.uplink,
								Path_Delay, Path_Average_Atmp, aglT, An, Bn, Wn, RndPhase, aglR,
								distypeR, distypeT,
								&para.UpdateInterval, &para.UpdatesPerBurst, &para.dt, &para.NumOfTaps, 
								&para.aglsprR, &para.aglsprT, &para.dR, &para.dT, 
								&para.corrmodel, &para.alf, &para.fd, &para.loson);



				int len_Tx = Data_Tx.size()*para.Nt;
				int len_Tx_2 = Data_Tx_2.size();
				int len_An = An.size();
				int len_Bn = Bn.size();
				int len_Wn = Wn.size();
				int len_RndPhase = RndPhase.size();
				int len_Delay = Path_Delay.size();
				int len_Average_Atmp = Path_Average_Atmp.size();
				int len_aglT = aglT.size();
				int len_aglR = aglR.size();
				int len_distypeT = distypeT.size();
				int len_distypeR = distypeR.size();
				int len_foreData = para.Nr*Path_Delay.back();

				double* ch_Data_Tx_pr = new double[len_Tx];
				double* ch_Data_Tx_pi = new double[len_Tx];
				//double* ch_Data_Tx_2_pr = new double[len_Tx_2];
				//double* ch_Data_Tx_2_pi = new double[len_Tx_2];
				double* ch_An = new double[len_An];
				double* ch_Bn = new double[len_Bn];
				double* ch_Wn = new double[len_Wn];
				double* ch_RndPhase = new double[len_RndPhase];
				double* ch_Delay = new double[len_Delay];
				double* ch_Average_Atmp = new double[len_Average_Atmp];
				double* ch_aglT = new double[len_aglT];
				double* ch_aglR = new double[len_aglR];
				double* ch_distypeT = new double[len_distypeT];
				double* ch_distypeR = new double[len_distypeR];
				double* ch_fore_data_pr = new double[len_foreData];
				double* ch_fore_data_pi = new double[len_foreData];

				double* ch_Data_Out_pr = new double[len_Tx];
				double* ch_Data_Out_pi = new double[len_Tx];
				for (int i=0; i<len_Tx; i++)
				{
					ch_Data_Out_pr[i]=0;
					ch_Data_Out_pi[i]=0;
				}

				double* ch_Fading_Out_pr = new double[len_Tx*para.NumOfTaps*para.Nt*para.Nr];
				double* ch_Fading_Out_pi = new double[len_Tx*para.NumOfTaps*para.Nt*para.Nr];

				switch(para.Nt)
				{
				case 1:
					for(int j=0;j<len_Tx;j++) ch_Data_Tx_pr[j]=Data_Tx[j].real();
					for(int j=0;j<len_Tx;j++) ch_Data_Tx_pi[j]=Data_Tx[j].imag();
					break;
				case 2:
					for(int j=0;j<(len_Tx/2);j++) ch_Data_Tx_pr[2*j]=Data_Tx[j].real();
					for(int j=0;j<(len_Tx/2);j++) ch_Data_Tx_pr[2*j+1]=Data_Tx_2[j].real();
					for(int j=0;j<(len_Tx/2);j++) ch_Data_Tx_pi[2*j]=Data_Tx[j].imag();
					for(int j=0;j<(len_Tx/2);j++) ch_Data_Tx_pi[2*j+1]=Data_Tx_2[j].imag();
					break;
				}

				for(int j=0;j<len_An;j++) ch_An[j]=An[j];
				for(int j=0;j<len_Bn;j++) ch_Bn[j]=Bn[j];
				for(int j=0;j<len_Wn;j++) ch_Wn[j]=Wn[j];
				for(int j=0;j<len_RndPhase;j++) ch_RndPhase[j]=RndPhase[j];
				for(int j=0;j<len_Delay;j++) ch_Delay[j]=Path_Delay[j];
				for(int j=0;j<len_Average_Atmp;j++) ch_Average_Atmp[j]=Path_Average_Atmp[j];
				for(int j=0;j<len_aglT;j++) ch_aglT[j]=aglT[j];
				for(int j=0;j<len_aglR;j++) ch_aglR[j]=aglR[j];
				for(int j=0;j<len_distypeT;j++) ch_distypeT[j]=distypeT[j];
				for(int j=0;j<len_distypeR;j++) ch_distypeR[j]=distypeR[j];
				for(int j=0;j<len_foreData;j++) ch_fore_data_pr[j]=0.0;
				for(int j=0;j<len_foreData;j++) ch_fore_data_pi[j]=0.0;

				ch_fore_data_pr[0] = 0.0000000001;
				ch_fore_data_pi[0] = 0.0000000001;

				////Data_Tx过信道
				corrj_channel(	ch_Data_Tx_pr, ch_Data_Tx_pi, ch_An, ch_Bn, ch_Wn, 
							ch_RndPhase, para.NumofFreq, para.UpdateInterval, para.UpdatesPerBurst, 
							para.LengthOfFrame, para.SampleIndex, para.dt, para.NumOfTaps, 
							ch_Delay, ch_Average_Atmp,  
							ch_fore_data_pr,ch_fore_data_pi, ch_Data_Out_pr, 
							ch_Data_Out_pi, ch_Fading_Out_pr, ch_Fading_Out_pi,
							para.Nr,para.Nt,para.aglsprR,ch_aglR,para.aglsprT,ch_aglT,
							para.dR,para.dT,ch_distypeR,ch_distypeT,para.corrmodel,para.alf,
							para.fd,para.loson,para.uplink);


				Data_Tx.clear();
				Data_Tx_2.clear();
                


				switch(para.Nr)
					{
					case 1:
						for(int j=0;j<len_Tx;j++)Data_Tx.push_back(complex<double>::complex(ch_Data_Out_pr[j],ch_Data_Out_pi[j]));
						break;
					case 2:
						for(int j=0;j<(len_Tx/2);j++) Data_Tx.push_back(complex<double>::complex(ch_Data_Out_pr[2*j],ch_Data_Out_pi[2*j]));
						for(int j=0;j<(len_Tx/2);j++) Data_Tx_2.push_back(complex<double>::complex(ch_Data_Out_pr[2*j+1],ch_Data_Out_pi[2*j+1]));
						break;
					}


				


				delete[] ch_Data_Tx_pr;
				delete[] ch_Data_Tx_pi;
				delete[] ch_Data_Out_pr;
				delete[] ch_Data_Out_pi;
				delete[] ch_Fading_Out_pr;
				delete[] ch_Fading_Out_pi;

				delete[] ch_An;
				delete[] ch_Bn;
				delete[] ch_Wn;
				delete[] ch_RndPhase;
				delete[] ch_Delay;
				delete[] ch_Average_Atmp;
				delete[] ch_aglT;
				delete[] ch_aglR;
				delete[] ch_distypeT;
				delete[] ch_distypeR;
				if (para.Channel_Type < 4)//模式4为单径，无需删除该数组
				{
					delete[] ch_fore_data_pr;
					delete[] ch_fore_data_pi;
				}
				}
			    /**********************************************
								 Oversampling
                ***********************************************/
				vector<double> filter_coef;
				vector<complex<double>> filter_coef_temp;
			//	double filter_alpha = 1.5/double(para.Oversampling_Ratio);
				double filter_alpha=0.24;
				int filter_length = 30;
				//if(para.Oversampling_Ratio > 1) 
					filter_coef = Root_raised_cosine(filter_alpha,(double)para.Oversampling_Ratio, filter_length);//生成窗函数
				
				//for (int i=0; i < filter_length ; i++)
				//{
				//	if (i<3||i>28)
				//	{
				//			filter_coef_temp.push_back(0.1);
				//	} 
				//	else
				//	{
				//			filter_coef_temp.push_back(1);
				//	}


				//}

				//	IFFT(filter_coef_temp, filter_length);

				//for (int i=0; i < filter_length ; i++)
				//{
				//	filter_coef.push_back(filter_coef_temp[i].real());
				//}
				
				OverSample_Zero(Data_Tx, para.Oversampling_Ratio);//过采样补零
				if(para.Oversampling_Ratio > 1) 
				{
					My_conv(Data_Tx, filter_coef, Data_Tx.size());
				}
				if(para.Nr == 2) 
				{
					OverSample_Zero(Data_Tx_2, para.Oversampling_Ratio);
					if(para.Oversampling_Ratio > 1) 
						My_conv(Data_Tx_2, filter_coef, Data_Tx_2.size());
				}
				filter_coef.clear();
            }

void ProtocolE::Preamble_generate(vector <complex <double> > &Preamble_out)
{
// 函数名: Preamble_generate
// 函数功能描述:产生前导，一个OFDMA符号
// 输入参数: 
// 输出参数:Preamble_out
// Called By: tx_param
// 补充说明:
// 修改日期: 2017/08/07：15:25:42
	vector <double> Preamble_data;
	vector <double> Preamble_dec;
	vector <int> Preamble_data0, Preamble_data1, Preamble_data2, Preamble_data3, Preamble_data4, Preamble_data5, Preamble_data6, Preamble_data7, Preamble_data8, Preamble_data9, Preamble_data10, Preamble_data11, Preamble_data12, Preamble_data13, Preamble_data14, Preamble_data15;
	switch (para.FFT_Size)
	{
	case 512:
		hex2bit("66C9CB4D1", Preamble_data0);
		hex2bit("C8F31D60F", Preamble_data1);
		hex2bit("5795886EE", Preamble_data2);
		hex2bit("02FFF6BE4", Preamble_data3);

		Preamble_data.insert(Preamble_data.end(), Preamble_data0.begin(), Preamble_data0.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data1.begin(), Preamble_data1.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data2.begin(), Preamble_data2.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data3.begin(), Preamble_data3.end());
		Preamble_data0.clear();
		Preamble_data1.clear();
		Preamble_data2.clear();
		Preamble_data3.clear();

		Preamble_out.resize(512);
		for(int i=0;i<143;i++)
		{
			Preamble_out[para.segment_ID+3*i+42]=complex<double>(4*(double)sqrt((double) 2)*(0.5-Preamble_data[i]),0);
		}
		Preamble_out[256]=complex<double>(0,0);
		Preamble_data.clear();
		break;
	case 1024:
		hex2bit("A6F294537", Preamble_data0);
		hex2bit("B285E1844", Preamble_data1);
		hex2bit("677D133E4", Preamble_data2);
		hex2bit("D53CCB1F1", Preamble_data3);
		hex2bit("82DE00489", Preamble_data4);
		hex2bit("E53E6B6E7", Preamble_data5);
		hex2bit("7065C7EE7", Preamble_data6);
		hex2bit("D0ADBEAF", Preamble_data7);

		Preamble_data.insert(Preamble_data.end(), Preamble_data0.begin(), Preamble_data0.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data1.begin(), Preamble_data1.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data2.begin(), Preamble_data2.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data3.begin(), Preamble_data3.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data4.begin(), Preamble_data4.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data5.begin(), Preamble_data5.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data6.begin(), Preamble_data6.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data7.begin(), Preamble_data7.end());
		Preamble_data0.clear();
		Preamble_data1.clear();
		Preamble_data2.clear();
		Preamble_data3.clear();
		Preamble_data4.clear();
		Preamble_data5.clear();
		Preamble_data6.clear();
		Preamble_data7.clear();

		Preamble_out.resize(1024);
		for(int i=0;i<284;i++)
		{
			Preamble_out[para.segment_ID+3*i+86]=complex<double>(4*(double)sqrt((double) 2)*(0.5-Preamble_data[i]),0);
		}
		Preamble_out[512]=complex<double>(0,0);

		Preamble_data.clear();
		break;
	case 2048:
		hex2bit("C12B7F736", Preamble_data0);
		hex2bit("CFFB14B6A", Preamble_data1);
		hex2bit("BF4EB50A6", Preamble_data2);
		hex2bit("0B7A3B416", Preamble_data3);
		hex2bit("3EA3360F6", Preamble_data4);
		hex2bit("97C450759", Preamble_data5);
		hex2bit("97ACE17BB", Preamble_data6);
		hex2bit("1512C7C0C", Preamble_data7);
		hex2bit("EBB34B389", Preamble_data8);
		hex2bit("D8784553C", Preamble_data9);
		hex2bit("0FC60BDE4", Preamble_data10);
		hex2bit("F166CF7B0", Preamble_data11);
		hex2bit("4856442D9", Preamble_data12);
		hex2bit("7539FB915", Preamble_data13);
		hex2bit("D80820CED", Preamble_data14);
		hex2bit("D858483", Preamble_data15);

		Preamble_data.insert(Preamble_data.end(), Preamble_data0.begin(), Preamble_data0.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data1.begin(), Preamble_data1.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data2.begin(), Preamble_data2.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data3.begin(), Preamble_data3.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data4.begin(), Preamble_data4.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data5.begin(), Preamble_data5.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data6.begin(), Preamble_data6.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data7.begin(), Preamble_data7.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data8.begin(), Preamble_data8.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data9.begin(), Preamble_data9.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data10.begin(), Preamble_data10.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data11.begin(), Preamble_data11.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data12.begin(), Preamble_data12.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data13.begin(), Preamble_data13.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data14.begin(), Preamble_data14.end());
		Preamble_data.insert(Preamble_data.end(), Preamble_data15.begin(), Preamble_data15.end());
		Preamble_data0.clear();
		Preamble_data1.clear();
		Preamble_data2.clear();
		Preamble_data3.clear();
		Preamble_data4.clear();
		Preamble_data5.clear();
		Preamble_data6.clear();
		Preamble_data7.clear();
		Preamble_data8.clear();
		Preamble_data9.clear();
		Preamble_data10.clear();
		Preamble_data11.clear();
		Preamble_data12.clear();
		Preamble_data13.clear();
		Preamble_data14.clear();
		Preamble_data15.clear();
		Preamble_out.resize(2048);
		for(int i=0;i<568;i++)
		{
			Preamble_out[para.segment_ID+3*i+172]=complex<double>(4*(double)sqrt((double) 2)*(0.5-Preamble_data[i]),0);
		}
		Preamble_out[1024]=complex<double>(0,0);
		Preamble_data.clear();
		break;
	}
}
void ProtocolE::SYNC_generate(vector <complex <double> > &SYNC_out)
{
	vector <int> SYNC_Modulated, SYNC_data;
	vector <double> SYNC_dec;
	switch (para.FFT_Size)
	{
	case 128:
		hex2bit("590A18B643F9D0", SYNC_data);

		SYNC_Modulated.resize(128);
		for(int i=1;i<=53;i++)
		{
			SYNC_Modulated[2*i-1+11]=1-2*SYNC_data[i];
		}

		SYNC_data.clear();
	case 512:
		hex2bit("5642862D90FE75642862A6F018B642862D90FE749BD79D590FE740", SYNC_data);

		SYNC_Modulated.resize(512);
		for(int i=1;i<=213;i++)
		{
			SYNC_Modulated[2*i-1+43]=1-2*SYNC_data[i];
		}

		SYNC_data.clear();
		break;
	case 1024:
		hex2bit("473A0B21CE9537F3A0B20316AC873A0B21CE95378C5F4DFCE9537F3A0B21CE9537F3A0B20316AC80C5F4DE316AC873A0B20316AC800", SYNC_data);

		SYNC_Modulated.resize(1024);
		for(int i=1;i<=425;i++)
		{
			SYNC_Modulated[2*i-1+87]=1-2*SYNC_data[i];
		}

		SYNC_data.clear();
		break;
		//case 2048:  //协议标准中，没有2048点的SYNC序列
		//	hex2bit("473A0B21CE9537F3A0B20316AC873A0B21CE95378C5F4DFCE9537F3A0B21CE9537F3A0B20316AC80C5F4DE316AC873A0B20316AC800", SYNC_data);

		//	SYNC_Modulated.resize(1024);
		//	for(int i=1;i<=425;i++)
		//	{
		//		SYNC_Modulated[2*i-1+87]=1-2*SYNC_data[i];
		//	}

		//	SYNC_data.clear();
		//break;
	}

	SYNC_out.insert(SYNC_out.end(),SYNC_Modulated.begin(),SYNC_Modulated.end());
}
void ProtocolE::poilts_generate(vector <int> &Poilts_Out)
{
// 函数名: poilts_generate
// 函数功能描述:产生BPSK调制后的导频序列
// 输入参数: 
// 输出参数:Poilts_Out
// Called By: tx_param
// 补充说明:对于下行PRBS_ID=00000 01 1111 b0-b4是IDCell，本文中是0，b5-b7=segment+1，segment=0,b7-b10=1111
// 修改日期: 2017/08/07：15:30:58
	vector <int> Poilts, Poilts_New, DataSubC_Sym;
	int feedback;
	Poilts.push_back(0);
	Poilts.push_back(0);
	Poilts.push_back(0);
	Poilts.push_back(0);
	Poilts.push_back(0);
	Poilts.push_back(0);
	Poilts.push_back(1);
	Poilts.push_back(1);
	Poilts.push_back(1);
	Poilts.push_back(1);
	Poilts.push_back(1);

	for (int i = 0; i < 3000; i++)
	{
		feedback = Poilts[8] ^ Poilts[10];
		Poilts_Out.push_back(1-2*Poilts[10]);
		circshift(Poilts, 1, Poilts_New);
		for (int j = 0; j < Poilts.size(); j++)
		{
			Poilts[j] = Poilts_New[j];
		}
		Poilts_New.clear();
		Poilts[0] = feedback;
	}
}
void ProtocolE::randomization_generate( int UncBSize,vector <int> &RadomS)
{
// 函数名: randomization_generate
// 函数功能描述:产生扰码序列
// 输入参数: UncBSize
// 输出参数:RadomS
// Called By: channel_code
// 补充说明:前导和FCH不进行扰码 其余均要进行扰码
// 修改日期: 2017/08/07：15:37:06
	vector <int> IniVec,IniVec_tmp;
	IniVec.push_back(0);
	IniVec.push_back(1);
	IniVec.push_back(1);
	IniVec.push_back(0);
	IniVec.push_back(1);
	IniVec.push_back(1);
	IniVec.push_back(1);
	IniVec.push_back(0);
	IniVec.push_back(0);
	IniVec.push_back(0);
	IniVec.push_back(1);
	IniVec.push_back(0);
	IniVec.push_back(1);
	IniVec.push_back(0);
	IniVec.push_back(1);

	for (int ii = 0; ii < UncBSize; ii++)
	{
		RadomS.push_back(IniVec[13] ^ IniVec[14]);
		circshift(IniVec, 1, IniVec_tmp);
		IniVec_tmp[0] = RadomS[ii];
		for (int i = 0; i < IniVec.size(); i++)
		{
			IniVec[i] = IniVec_tmp[i];
		}
		IniVec_tmp.clear();
	}
}
void ProtocolE::interleaving_generate(vector <int> &data_in,int Modulation_Mode,int Ncbpss,vector <int> &InterleavedData)
{
// 函数名: interleaving_generate
// 函数功能描述:对数据进行交织
// 输入参数: data_in、Modulation_Mode、Ncbpss
// 输出参数:InterleavedData
// Called By: channel_code、FCH_modulation
// 补充说明:
// 修改日期: 2017/08/07：15:39:50
	int s;
	s = ceil(((float)log((float)Modulation_Mode) / (float)log((float)2))/2.0); //ceil正方向舍入

	vector<int> firstStepIndex(Ncbpss);
	vector<int> secondStepIndex(Ncbpss);
	vector<int> k;
	for(int i=0;i<Ncbpss;++i)
		k.push_back(i);
	int* dataTemp1 = new int [Ncbpss];
	int* dataInterleaved = new int [Ncbpss];
	for(int i=0;i<Ncbpss;++i)
	{
		firstStepIndex[i] = (Ncbpss/16)*(k[i]%16) + k[i]/16;
	//	dataTemp1[firstStepIndex[i]] = data_in[i];
	}
	//Second permutation
	for(int i=0;i<Ncbpss;++i)
	{
		secondStepIndex[i] = s*(firstStepIndex[i]/s) + (firstStepIndex[i]+Ncbpss-16*firstStepIndex[i]/Ncbpss)%s;
		dataInterleaved[secondStepIndex[i]] = data_in[i];
	}
	for(int i=0;i<Ncbpss;i++)
		InterleavedData.push_back(dataInterleaved[i]);

}

void ProtocolE::FCH_generate(vector <int> &DL_Frame_Prefix ,int* Group)
{
// 函数名: FCH_generate
// 函数功能描述:产生DL_Frame_Prefix，DL_Frame_Prefix重复4倍即为FCH
// 输入参数: Group、全局参数para.DL_MAP_Length、para.DL_MAP_Coding、para.DL_MAP_repetition
// 输出参数:DL_Frame_Prefix
// Called By: tx_param
// 补充说明:
// 修改日期: 2017/08/07：15:41:26
		vector <int> a_binary;;
		for (int i = 0; i < 24; i++)
		{
			DL_Frame_Prefix.push_back(0);
		}
		if (Group[0] = 0)
			DL_Frame_Prefix[0] = 0;
		else
			DL_Frame_Prefix[0] = 1;
		if (Group[1] = 0)
			DL_Frame_Prefix[1] = 0;
		else
			DL_Frame_Prefix[1] = 1;
		if (Group[2] = 0)
			DL_Frame_Prefix[2] = 0;
		else
			DL_Frame_Prefix[2] = 1;
		if (Group[3] = 0)
			DL_Frame_Prefix[3] = 0;
		else
			DL_Frame_Prefix[3] = 1;
		if (Group[4] = 0)
			DL_Frame_Prefix[4] = 0;
		else
			DL_Frame_Prefix[4] = 1;
		if (Group[5] = 0)
			DL_Frame_Prefix[5] = 0;
		else
			DL_Frame_Prefix[5] = 1;

		switch (para.DL_MAP_repetition)
		{
		case 0:
			DL_Frame_Prefix[7] = 0;
			DL_Frame_Prefix[8] = 0;
			break;
		case 2:
			DL_Frame_Prefix[7] = 0;
			DL_Frame_Prefix[8] = 1;
			break;
		case 4:
			DL_Frame_Prefix[7] = 1;
			DL_Frame_Prefix[8] = 0;
			break;
		case 6:
			DL_Frame_Prefix[7] = 1;
			DL_Frame_Prefix[8] = 1;
			break;
		}

		switch (para.DL_MAP_Coding)
		{
		case 0:
			DL_Frame_Prefix[9] = 0;
			DL_Frame_Prefix[10] = 0;
			DL_Frame_Prefix[11] = 0;
			break;
		case 1:
			DL_Frame_Prefix[9] = 0;
			DL_Frame_Prefix[10] = 0;
			DL_Frame_Prefix[11] = 1;
			break;
		case 2:
			DL_Frame_Prefix[9] = 0;
			DL_Frame_Prefix[10] = 1;
			DL_Frame_Prefix[11] = 0;
			break;
		case 3:
			DL_Frame_Prefix[9] = 0;
			DL_Frame_Prefix[10] = 1;
			DL_Frame_Prefix[11] = 1;
			break;
		case 4:
			DL_Frame_Prefix[9] = 1;
			DL_Frame_Prefix[10] = 0;
			DL_Frame_Prefix[11] = 0;
			break;
		}
		de2bi(para.DL_MAP_Length, 8, a_binary);
		for (int i = 0; i < 8; i++)
		{
			DL_Frame_Prefix[19 - i] = a_binary[i];
		}

}
int  ProtocolE::deg(int poly[],int t)
{
// 函数名:deg 
// 函数功能描述:计算多项式次数，用来产生HCS
// 输入参数: poly、t
// 输出参数:返回p
// Called By: FCH_generate
// 补充说明:
// 修改日期: 2017/08/07：14:52:34
	int l,p=0;
	for(l=0;l<t;l++)
	{
		if(poly[l])
		{
			if(l>=p) p=l;
		}
	}
	return p;
}
void ProtocolE::MAC_head(int CI,int length,int CID,vector <int> &MAC_head_DLMAP)
{
// 函数名: MAC_head
// 函数功能描述:产生MAC帧头，与DL_MAP共同构成DL_MAP PDU
// 输入参数: CI、length、CID
// 输出参数:MAC_head_DLMAP
// Called By: tx_param
// 补充说明:DL_MAP模式为normal，且无CRC校验
// 修改日期: 2017/08/07：15:44:07
	int mac_head[48]={0};
	vector <int> a_binary;

	de2bi(CI, 1, a_binary);
	mac_head[9]=a_binary[0];
	a_binary.clear();

	de2bi(length, 11, a_binary);
	for (int i = 0; i < 11; i++)
	{
		mac_head[23 - i] = a_binary[i];
	}
	a_binary.clear();

	de2bi(CID, 16, a_binary);
	for (int i = 0; i < 16; i++)
	{
		mac_head[39 - i] = a_binary[i];
	}
	int N=100;
	int a[100]={0},b[100]={0},bcs[100]={0},cs[100]={0},q[100]={0};
	int temp;
	for(int i=0;i<40;i++)
		a[8+39-i]=mac_head[i];
	b[8]=b[2]=b[1]=b[0]=1;
	
	for(int i=0;i<N;i++)
	{
		bcs[i]=a[i];
		cs[i]=b[i];
	}
	temp=deg(bcs,N)-deg(cs,N);
	while(temp>=0)
	{
		q[temp]=1;
		for(int j=0;j<=deg(cs,N);j++)
			bcs[temp+j]^=b[j];
		temp=deg(bcs,N)-deg(cs,N);
	}

	for(int i=0;i<8;i++)
	{
		if(bcs[i]==1) 
		mac_head[47-i]=bcs[i];
	}
	MAC_head_DLMAP.insert(MAC_head_DLMAP.end(), mac_head, mac_head + 48);
}

void ProtocolE::DL_MAP_generate(int Power_boost, int DIUC_index, int OFDMA_Symbol_offset, int Subchannel_offset0, int Subchannel_offset, int Num_OFDMA_Symbols, int Usedsubchannel, int RCI, vector <int> &DL_MAP)
{
// 函数名: DL_MAP_generate
// 函数功能描述:产生下行突发的映射信息
// 输入参数: Power_boost、DIUC_index、OFDMA_Symbol_offset、Subchannel_offset、Num_OFDMA_Symbols、Usedsubchannel、RCI、
// 输出参数:DL_MAP
// Called By: tx_param
// 补充说明:DL_MAP中有两个DL_MAP_IE，第一个IE参数固定，第二个IE参数由形参决定
// 修改日期: 2017/08/07：15:51:15
	int Power_boost1[3];
	switch (Power_boost)
	{
	case 0:
		Power_boost1[0] = 0;
		Power_boost1[1] = 0;
		Power_boost1[2] = 0;
		break;
	case 6:
		Power_boost1[0] = 0;
		Power_boost1[1] = 0;
		Power_boost1[2] = 1;
		break;
	case -6:
		Power_boost1[0] = 0;
		Power_boost1[1] = 1;
		Power_boost1[2] = 0;
		break;
	case 9:
		Power_boost1[0] = 0;
		Power_boost1[1] = 1;
		Power_boost1[2] = 1;
		break;
	case 3:
		Power_boost1[0] = 1;
		Power_boost1[1] = 0;
		Power_boost1[2] = 0;
		break;
	case -3:
		Power_boost1[0] = 1;
		Power_boost1[1] = 0;
		Power_boost1[2] = 1;
		break;
	case -9:
		Power_boost1[0] = 1;
		Power_boost1[1] = 1;
		Power_boost1[2] = 0;
		break;
	case -12:
		Power_boost1[0] = 1;
		Power_boost1[1] = 1;
		Power_boost1[2] = 1;
		break;
	}
	int DL_MAP_IE[36] = {0};
	int DL_MAP_IE0[36] = {0};
	int DL_MAP_RE[104] = {0};
	int STC_DL_Zone_IE[44]={0};
	int DL_MAP_IE2[36]={0};
	vector <int> a_binary;
	de2bi(DIUC_index, 4, a_binary);
	for (int i = 0; i < 4; i++)
	{
		DL_MAP_IE[3 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(OFDMA_Symbol_offset, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_IE[11 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(Subchannel_offset, 6, a_binary);
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE[17 - i] = a_binary[i];
	}
	a_binary.clear();
	for (int i = 0; i < 3; i++)
	{
		DL_MAP_IE[18 + i] = Power_boost1[i];
	}
	de2bi(Num_OFDMA_Symbols, 7, a_binary);
	for (int i = 0; i < 7; i++)
	{
		DL_MAP_IE[27 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(Usedsubchannel, 6, a_binary);
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE[33 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(RCI/2, 2, a_binary);
	for (int i = 0; i < 2; i++)
	{
		DL_MAP_IE[35 - i] = a_binary[i];
	}
	a_binary.clear();
	//IE0
	de2bi(0, 4, a_binary); //DIUC
	for (int i = 0; i < 4; i++)
	{
		DL_MAP_IE0[3 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(1, 8, a_binary);  //符号偏移
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_IE0[11 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(Subchannel_offset0, 6, a_binary);//子信道偏移
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE0[17 - i] = a_binary[i];
	}
	a_binary.clear();
	for (int i = 0; i < 3; i++)
	{
		DL_MAP_IE0[18 + i] = Power_boost1[i];
	}
	de2bi(2, 7, a_binary);//突发所占符号数
	for (int i = 0; i < 7; i++)
	{
		DL_MAP_IE0[27 - i] = a_binary[i];
	}
	a_binary.clear();

	de2bi(para.Usedsubchannel-Subchannel_offset0, 6, a_binary);//子信道数

	
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE0[33 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(0, 2, a_binary);
	for (int i = 0; i < 2; i++)
	{
		DL_MAP_IE0[35 - i] = a_binary[i];
	}
	a_binary.clear();

	//STC_DL_Zone_IE
	de2bi(15, 4, a_binary);//DIUC
	for (int i = 0; i < 4; i++)
	{
		STC_DL_Zone_IE[3 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(1, 4, a_binary);//Extended DIUC
	for (int i = 0; i < 4; i++)
	{
		STC_DL_Zone_IE[7 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(4, 4, a_binary);//length
	for (int i = 0; i < 4; i++)
	{
		STC_DL_Zone_IE[11 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(5, 8, a_binary);//符号偏移
	for (int i = 0; i < 8; i++)
	{
		STC_DL_Zone_IE[19 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(0, 2, a_binary);//permutation
	for (int i = 0; i < 2; i++)
	{
		STC_DL_Zone_IE[21 - i] = a_binary[i];
	}
	a_binary.clear();
	STC_DL_Zone_IE[22]=1;//use all sc
	de2bi(1, 2, a_binary);//STC
	for (int i = 0; i < 2; i++)
	{
		STC_DL_Zone_IE[24 - i] = a_binary[i];
	}
	a_binary.clear();
	//IE1
	de2bi(0, 4, a_binary);
	for (int i = 0; i < 4; i++)
	{
		DL_MAP_IE2[3 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(5, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_IE2[11 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(0, 6, a_binary);
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE2[17 - i] = a_binary[i];
	}
	a_binary.clear();
	for (int i = 0; i < 3; i++)
	{
		DL_MAP_IE2[18 + i] = Power_boost1[i];
	}
	de2bi(10, 7, a_binary);
	for (int i = 0; i < 7; i++)
	{
		DL_MAP_IE2[27 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.Usedsubchannel, 6, a_binary);
	for (int i = 0; i < 6; i++)
	{
		DL_MAP_IE2[33 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(0, 2, a_binary);
	for (int i = 0; i < 2; i++)
	{
		DL_MAP_IE2[35 - i] = a_binary[i];
	}
	a_binary.clear();

	//DL_MAP  104比特
	de2bi(para.Management_Message_Type, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_RE[7 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.Frame_Duration_Code, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_RE[15 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.Frame_Number_index, 24, a_binary);
	for (int i = 0; i < 24; i++)
	{
		DL_MAP_RE[39 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.DCD_Count_index, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_RE[47 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.BS_ID_index, 48, a_binary);
	for (int i = 0; i < 48; i++)
	{
		DL_MAP_RE[95 - i] = a_binary[i];
	}
	a_binary.clear();
	de2bi(para.NumOfDL_Frame_Symbol, 8, a_binary);
	for (int i = 0; i < 8; i++)
	{
		DL_MAP_RE[103 - i] = a_binary[i];
	}
	a_binary.clear();

	DL_MAP.insert(DL_MAP.end(), DL_MAP_RE, DL_MAP_RE + 104);
	DL_MAP.insert(DL_MAP.end(), DL_MAP_IE0, DL_MAP_IE0 + 36);
	DL_MAP.insert(DL_MAP.end(), DL_MAP_IE, DL_MAP_IE + 36);
	DL_MAP.insert(DL_MAP.end(), STC_DL_Zone_IE, STC_DL_Zone_IE + 44);
	DL_MAP.insert(DL_MAP.end(), DL_MAP_IE2, DL_MAP_IE2 + 36);
	for(int i=0;i<32;i++)
	{
		DL_MAP.push_back(1);
	}
}

void ProtocolE::FCH_modulation(vector <int> &DL_Frame_Prefix, vector <complex <double> > &data_repetition)
{
// 函数名: FCH_modulation
// 函数功能描述:对DL_Frame_Prefix重复四次后进行QPSK CC/CTC 1/2调制
// 输入参数: DL_Frame_Prefix
// 输出参数:data_repetition
// Called By: tx_param
// 补充说明:FCH不扰码，且FEC块为1个时隙
// 修改日期: 2017/08/07：15:55:20
	vector <int>RadomS,jkValue,CC_in;

	vector <int> Re_data_after_coding, data_CC_coding,InterleavedData;
	vector <complex <double> >data_after_modulation ;
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<24;j++)
		{
			Re_data_after_coding.push_back(DL_Frame_Prefix[j]);
		}
	}

	CC_coding(Re_data_after_coding, 0.5, data_CC_coding);
	
	interleaving_generate(data_CC_coding,4,96,InterleavedData);

	modulation(InterleavedData,0,data_after_modulation);

	for (int i = 0; i < 4; i++)
	{
		for(int j=0;j<48;j++)
			data_repetition.push_back(data_after_modulation[j]);
	}	
}
void ProtocolE::DL_MAP_modulation(vector <int> &DL_MAP, int DL_MAP_repetition, int DL_MAP_Coding, vector <complex <double> > &DL_MAP_out)
{
// 函数名: DL_MAP_modulation
// 函数功能描述:下行映射调制
// 输入参数: DL_MAP、DL_MAP_repetition、DL_MAP_Coding
// 输出参数:DL_MAP_out
// Called By: tx_param
// 补充说明:DL_MAP信息需要扰码，且如果CC编码，FEC块为5个时隙
// 修改日期: 2017/08/07：15:56:59
	vector <complex <double> > DL_MAP_mod;
	vector <int> DL_MAP_channelcode;
	int lengthofDL_MAP=DL_MAP.size();
	int mcs;
	if(DL_MAP_Coding==0) mcs=0;
	else  mcs=7;
	channel_code(DL_MAP,mcs,lengthofDL_MAP,DL_MAP_channelcode);
	modulation(DL_MAP_channelcode,0, DL_MAP_mod);
	repetition(DL_MAP_repetition,DL_MAP_mod,DL_MAP_out);
}

void ProtocolE::block_channel_code(int UncBSize, vector <int>&burst2,int j_block,int MCS_index,vector <int>&InterLeavData2)
{
// 函数名: block_channel_code
// 函数功能描述:对待信道编码的数据分成数据块（根据CC编码的时隙级联规则）
// 输入参数: UncBSize，burst2，j_block，MCS_index
// 输出参数:InterLeavData2
// Called By: run
// 补充说明:
// 修改日期: 2017/08/07：16:00:18
	int n_block=burst2.size()/UncBSize;
	int k_block=n_block/j_block;
	int m_block=n_block%j_block;
	vector<int>InterLeavData2_1,InterLeavData2_2,InterLeavData2_3;
	vector<int>burst2_1,burst2_2,burst2_3;
	if(m_block==0)
	{
		channel_code(burst2, MCS_index, j_block*UncBSize,InterLeavData2);
	}
	else
	{
		burst2_1.insert(burst2_1.end(),burst2.begin(),burst2.begin()+(k_block-1)*j_block*UncBSize);
		channel_code(burst2_1, MCS_index, j_block*UncBSize,InterLeavData2_1);

		burst2.erase(burst2.begin(),burst2.begin()+(k_block-1)*j_block*UncBSize);

		burst2_2.insert(burst2_2.end(),burst2.begin(),burst2.begin()+ceil((float)(j_block+m_block)/2)*UncBSize);
		channel_code(burst2_2, MCS_index,ceil((float)(j_block+m_block)/2)*UncBSize,InterLeavData2_2);

		burst2_3.insert(burst2_3.end(),burst2.end()-(j_block+m_block)/2*UncBSize,burst2.end());
		channel_code(burst2_3, MCS_index, (j_block+m_block)/2*UncBSize,InterLeavData2_3);

		InterLeavData2.insert(InterLeavData2.end(),InterLeavData2_1.begin(),InterLeavData2_1.end());
		InterLeavData2.insert(InterLeavData2.end(),InterLeavData2_2.begin(),InterLeavData2_2.end());
		InterLeavData2.insert(InterLeavData2.end(),InterLeavData2_3.begin(),InterLeavData2_3.end());
	}
}
void ProtocolE::channel_code(vector <int> &Data_In, int MCS_index, int UncBsize,vector <int> &Data_Out)
{
// 函数名: channel_code
// 函数功能描述:对分块的数据进行扰码、卷积编码、扰码操作
// 输入参数: Data_In、MCS_index、UncBsize
// 输出参数:Data_Out
// Called By: block_channel_code
// 补充说明:
// 修改日期: 2017/08/07：16:02:25
	int Modulation_Mode;
	float CodeRate;
	switch(MCS_index){
	case 0:
		CodeRate = 0.5;Modulation_Mode=4;break;
	case 1:
		CodeRate = 0.75;Modulation_Mode=4;break;
	case 2:
		CodeRate = 0.5;Modulation_Mode=16;break;
	case 3:
		CodeRate = 0.75;Modulation_Mode=16;break;
	case 4:
		CodeRate = 0.5;Modulation_Mode=64;break;
	case 5:
		CodeRate = 2/3.0;Modulation_Mode=64;break;
	case 6:
		CodeRate = 0.75;Modulation_Mode=64;break;
	case 7:
		CodeRate = 0.5;Modulation_Mode=4;break;
	case 8:
		CodeRate = 0.75;Modulation_Mode=4;break;
	case 9:
		CodeRate = 0.5;Modulation_Mode=16;break;
	case 10:
		CodeRate = 0.75;Modulation_Mode=16;break;
	case 11:
		CodeRate = 0.5;Modulation_Mode=64;break;
	case 12:
		CodeRate = 2/3.0;Modulation_Mode=64;break;
	}

	vector <int> RadomS;
	randomization_generate(UncBsize,RadomS);

	int Ncbpss=UncBsize/CodeRate;

	vector <int> FEC_DataIn, Data_Rd,CC_in,CC_DataEncode;
	
	int k_index;
	FEC_DataIn.insert(FEC_DataIn.end(),Data_In.begin(),Data_In.end());
	k_index = FEC_DataIn.size() / UncBsize;  //UncBsize*k_index

	int len=UncBsize;
	int len_coded=int((float)len/CodeRate);

	vector<int> InterleavedData_in,InterleavedData;
	

	for (int i = 0; i < k_index; i++)
	{
		int* CTC_in = new int[len];
		int* CTC_out = new int[len_coded];
		
		for (int j = 0; j < UncBsize; j++)
		{
			CC_in.push_back(FEC_DataIn[i*UncBsize+j] ^ RadomS[j]);
		}

		if(MCS_index<7)
		{        
			CC_coding(CC_in, CodeRate, CC_DataEncode);
			interleaving_generate(CC_DataEncode,Modulation_Mode,Ncbpss,Data_Out);
			
			CC_in.clear();
			CC_DataEncode.clear();
			
		}
		else
		{   
			for(int j=0;j<len;j++)
			CTC_in[j]=CC_in[j];

			ctc_encoder(CTC_in, CodeRate,len, CTC_out);

			for(int j=0;j<len_coded;j++)
				Data_Out.push_back(CTC_out[j]);

			delete[] CTC_in;
			delete[] CTC_out;
		} 
	}

}
void ProtocolE::CC_coding(vector <int> &data_af_random, float CC_rate, vector <int> &data_after_coding)
{
// 函数名: CC_coding
// 函数功能描述:卷积编码器对每一个FEC块进行编码，且编码率为1/2，约束长度K=7，对于X生成多项式G1=171(OCT),对于Y生成多项式G2=133(OCT),请选择编码率对CC编码进行删余:1：1/2  2: 2/3  3: 3/4
// 输入参数: data_af_random、CC_rate
// 输出参数:data_after_coding
// Called By: channel_code、FCH_modulation
// 补充说明:
// 修改日期: 2017/08/07：15:12:19
    vector <int> data_after_coding1;

    int len1, len2, CC_index;
    int len = data_af_random.size();
    int d[6] = {0};
    int d_temp[6] = {0};
    for (int i = 0; i < 6; i++)
    {
        d[i] = data_af_random[len - 1 - i]; //取随机序列后6个作为卷积编码器的初始状态
    }
    for (int i = 0; i < len; i++)
    {
        data_after_coding1.push_back((data_af_random[i] + d[0] + d[1] + d[2] + d[5]) % 2); //模2作异或运算X output
        data_after_coding1.push_back((data_af_random[i] + d[1] + d[2] + d[4] + d[5]) % 2); //Y output

        d_temp[0] = data_af_random[i];
        for (int i = 0; i < 5; i++)
        {
            d_temp[i + 1] = d[i]; //寄存器移位
        }
        for (int i = 0; i < 6; i++)
        {
            d[i] = d_temp[i]; //寄存器状态更新
        }
    }

	if (CC_rate == 0.5)
        CC_index = 1;
    else if (CC_rate-0.6666666<0.000001)
        CC_index = 2;
    else if (CC_rate == 0.75)
        CC_index = 3;
    else if (CC_rate == (float)5/6)
        CC_index = 4;


    switch (CC_index)
    {
    case 1:
        data_after_coding.insert(data_after_coding.end(), data_after_coding1.begin(), data_after_coding1.end());
        len1 = data_af_random.size(); //随机化后序列长度
        len2 = data_after_coding.size(); //卷积编码后序列长度
        break;
    case 2: //删余操作，每四个比特删除第三个比特1/2除以3/4
        for (int i = 0; i < 2 * len / 4; i++)
        {
            data_after_coding.push_back(data_after_coding1[4 * i]);
            data_after_coding.push_back(data_after_coding1[4 * i + 1]);
            data_after_coding.push_back(data_after_coding1[4 * i + 3]);
        }
        len2 = data_after_coding.size();
        len1 = len2 / 3 * 2;
        break;
    case 3:
        for (int i = 0; i < 2 * len / 6; i++)
        {
            data_after_coding.push_back(data_after_coding1[6 * i]);
            data_after_coding.push_back(data_after_coding1[6 * i + 1]);
            data_after_coding.push_back(data_after_coding1[6 * i + 3]);
            data_after_coding.push_back(data_after_coding1[6 * i + 4]);
        }
        len2 = data_after_coding.size();
        len1 = len2 / 4 * 3;
        break;
    case 4:
//        for (int i = 0; i < 2 * len / 10; i++)
//        {
//            data_after_coding.push_back(data_after_coding1[10 * i]);
//            data_after_coding.push_back(data_after_coding1[10 * i + 1]);
//        }
        break;
    }
}
void ProtocolE::ctc_encoder(int*Data_In, double code_rate, int len, int* Data_Out)
{
	//input format must be 1 
	register int i;
	int 	ii, jj, NN;//kk 原有
	int     fec_code_type;
	int 	tmp1 = 0, tmp2 = 0, shifter, parity;
	char	*qq;
	int bset0_len = len;
	int bset1_len = len;
	int perm_len = len;
	int out_len;
	int* bset0 = new int[bset0_len];//for Y1&Y2
	int* bset1 = new int[bset1_len];//for Y1&Y2
	int* perm_data = new int[perm_len];//perm_data for ctc_interleaver_2nd
	char    state;

	switch (len)
	{
	case 48:  {fec_code_type = 0; break; }
	case 72:  {fec_code_type = 1; break; }
	case 96:  {fec_code_type = 2; break; }
	case 144: {fec_code_type = 3; break; }
	case 192: {fec_code_type = 4; break; }
	case 216: {fec_code_type = 5; break; }
	case 240: {fec_code_type = 6; break; }
	case 288: {fec_code_type = 7; break; }
	case 360: {fec_code_type = 8; break; }
	case 384: {fec_code_type = 9; break; }
	case 432: {fec_code_type = 10; break; }
	case 480:{ fec_code_type = 11; break; }
	/*default:{
		cout << "data length is wrong" << endl;
		Sleep(3000);
		exit(0); }*/
	}


	//                            N/Nsub P0 Puncture    备忘：第一列为NN，第二至五列为P
	char  ctc_parameters[12][5] = { { 24, 5, 0, 0, 0 },
	{ 36, 11, 18, 0, 18 },
	{ 48, 13, 24, 0, 24 },
	{ 72, 11, 6, 0, 6 },
	{ 96, 7, 48, 24, 72 },
	{ 108, 11, 54, 56, 2 },
	{ 120, 13, 60, 0, 60 },
	{ 144, 17, 74, 72, 2 },
	{ 180, 11, 90, 0, 90 },
	{ 192, 11, 96, 48, 144 },
	{ 216, 13, 108, 0, 108 },
	{ 240, 13, 120, 60, 180 } };

	char Npsub[7] = { 24, 24, 24, 48, 48, 72, 72 };
	//This table is worked well with the current cmodel.
	//Nov/22/05 Miffie
	/*char sc[7][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },
	{ 0, 3, 1, 2, 7, 4, 6, 5 },
	{ 0, 2, 6, 4, 5, 7, 3, 1 },
	{ 0, 5, 7, 2, 6, 3, 1, 4 },
	{ 0, 6, 3, 5, 7, 1, 4, 2 },
	{ 0, 7, 1, 6, 3, 4, 2, 5 },
	{ 0, 3, 4, 7, 1, 2, 5, 6 } };*/
		char sc[7][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },
	{ 0, 6, 4, 2, 7, 1, 3, 5 },
	{ 0, 3, 7, 4, 5, 6, 2, 1 },
	{ 0, 5, 3, 6, 2, 7, 1, 4 },
	{ 0, 4, 1, 5, 6, 2, 7, 3 },
	{ 0, 2, 5, 7, 1, 3, 4, 6 },
	{ 0, 7, 6, 1, 3, 4, 5, 2 } };

	int si_parameters[12][2] = { { 3, 3 },//N=24
	{ 4, 3 },//N=36
	{ 4, 3 },//N=48
	{ 5, 3 },//N=72
	{ 5, 3 },//N=96
	{ 5, 4 },//N=108
	{ 6, 2 },//N=120
	{ 6, 3 },//N=144
	{ 6, 3 },//N=180
	{ 6, 3 },//N=192
	{ 6, 4 },//N=216
	{ 7, 2 } };//N=240 


	NN = ctc_parameters[fec_code_type][0];

	
	for (i = 0; i < len; i++)  bset0[i] = Data_In[i];
	for (i = 0; i < len; i++)  bset1[i] = Data_In[i];
	ctc_interleaver(bset1, bset1_len, ctc_parameters[fec_code_type][1],
		ctc_parameters[fec_code_type][2], ctc_parameters[fec_code_type][3],
		ctc_parameters[fec_code_type][4], bset1);
	
	//C1 Constituent Encoder
	state = 0;
	for (ii = 0; ii<bset0_len; ii += 2) { //C1
		constituent_encoder(bset0[ii], bset0[ii + 1],
			  &state, &tmp1, &tmp2);
	} //C1
	//cout << "C1 last state=" << endl;
	//cout << state;
	state = sc[NN % 7][state];
	for (ii = 0; ii<bset0_len; ii += 2) { //C1
		constituent_encoder(bset0[ii], bset0[ii + 1],
			&state, &bset0[ii], &bset0[ii + 1]);
	} //C1
	

	//C2 Constituent Encoder
	state = 0;
	for (ii = 0; ii<bset1_len; ii += 2) { //C2
		constituent_encoder(bset1[ii], bset1[ii + 1],
			&state, &tmp1, &tmp2);
	} //C2
	//cout << "C2 last state=" << endl;
	//cout << state << endl;
	state = sc[NN % 7][state];
	for (ii = 0; ii<bset1_len; ii += 2) { //C2
		constituent_encoder(bset1[ii], bset1[ii + 1],
			 &state, &bset1[ii], &bset1[ii + 1]);
	} //C2

	vector <int> A,B,Y1,Y2,W1,W2;
	vector <int> AB,Y1Y2,W1W2;

	//2nd interleaving
	//data 2nd interleaving
	memset(perm_data, 0, sizeof(int)*len);//清零
	subblock_interleaver(Data_In, len, 
		si_parameters[fec_code_type][0], si_parameters[fec_code_type][1], perm_data);
	for (i = 0; i < len/2; i++)  
	{
		A.push_back(perm_data[i*2]) ;
		B.push_back(perm_data[i*2+1]) ;
	}
		
	memset(perm_data, 0, sizeof(int)*len);

	//Y1Y2 2nd interleaving
	subblock_interleaver(bset0, bset0_len,
		si_parameters[fec_code_type][0], si_parameters[fec_code_type][1], perm_data);
	for (i = 0; i < len/2; i++)  
	{
		Y1.push_back(perm_data[i*2]) ;
		W1.push_back(perm_data[i*2+1]) ;
	}

//	symbol_grouping(perm_data, bset0_len, bset0);
	memset(perm_data, 0, sizeof(int)*bset0_len);

	//W1W2 2nd interleaving
	subblock_interleaver(bset1, bset1_len,
		si_parameters[fec_code_type][0], si_parameters[fec_code_type][1], perm_data);
	for (i = 0; i < len/2; i++)  
	{
		Y2.push_back(perm_data[i*2]) ;
		W2.push_back(perm_data[i*2+1]) ;
	}

	memset(perm_data, 0, sizeof(int)*bset1_len);



//bit grouping
	for (i = 0; i < len/2; i++) 
		AB.push_back(A[i]);
	for (i = 0; i < len/2; i++) 
		AB.push_back(B[i]);
	for (i = 0; i < len/2; i++)
	{
		Y1Y2.push_back(Y1[i]);
		Y1Y2.push_back(Y2[i]);
	}
	for (i = 0; i < len/2; i++)
	{
		W1W2.push_back(W1[i]);
		W1W2.push_back(W2[i]);
	}

	//Puncturing
	if (code_rate == 0.5) { 
		out_len = len * 2;
		for (ii = 0; ii<len; ii += 1) { 
			Data_Out[ii] = AB[ii];
		} 
		for (ii = 0; ii<(out_len - len); ii += 1) { 
			Data_Out[ii + len] = Y1Y2[ii];
		} 
	} 
	else if (code_rate-0.666666<0.0001) { 
		out_len = len * 3 / 2;
		for (ii = 0; ii<len; ii += 1) { 
			Data_Out[ii] = AB[ii];
		} 
		for (ii = 0; ii<(out_len - len); ii += 1) { 
			Data_Out[ii + len] = Y1Y2[ii];
		} 
	} 
	else if (code_rate == .75) { 
		out_len = len * 4 / 3;
		for (ii = 0; ii<len; ii += 1) { 
			Data_Out[ii] = AB[ii];
		} 
		for (ii = 0; ii<(out_len - len); ii += 1) { 
			Data_Out[ii + len] = Y1Y2[ii];
		} 
	} 
	else if (code_rate == (float)5/6) { //4:5/6
		out_len = len * 6 / 5;
		for (ii = 0; ii<len; ii += 1) { //4:5/6
			Data_Out[ii] = AB[ii];
		} //4:5/6
		for (ii = 0; ii<(out_len - len); ii += 1) { //4:5/6
			Data_Out[ii + len] = Y1Y2[ii];
		} //4:5/6
	} //4:5/6
	delete[] bset0;
	delete[] bset1;
	delete[] perm_data;
} //ctc_encoder
void ProtocolE::modulation (vector <int> &data_for_mod, int mcs_index,vector <complex <double> > &mod_out)
{
// 函数名: modulation
// 函数功能描述:对编码数据进行调制
// 输入参数: data_for_mod、mcs_index
// 输出参数:mod_out
// Called By: run、FCH_modulation
// 补充说明:
// 修改日期: 2017/08/07：14:30:58
	int modulation_mode;
	switch(mcs_index)
	{
	case 0:
		modulation_mode = 4;break;
	case 1:
		modulation_mode = 4;break;
	case 2:
		modulation_mode = 16;break;
	case 3:
		modulation_mode = 16;break;
	case 4:
		modulation_mode = 64;break;
	case 5:
		modulation_mode = 64;break;
	case 6:
		modulation_mode = 64;break;
	case 7:
		modulation_mode = 4;break;
	case 8:
		modulation_mode = 4;break;
	case 9:
		modulation_mode = 16;break;
	case 10:
		modulation_mode = 16;break;
	case 11:
		modulation_mode = 64;break;
	case 12:
		modulation_mode = 64;break;
	}
	switch(modulation_mode)
	{
	case 4:

		for (int i = 0; i < data_for_mod.size() / 2; i++)
		{
			mod_out.push_back(QPSK[2*data_for_mod[i*2]+data_for_mod[i*2+1]]);
		}
		break;
	case 16:
		for(int i=0;i<data_for_mod.size()/4;i++)
		{
			mod_out.push_back(_16QAM[8*data_for_mod[i*4]+4*data_for_mod[i*4+1]+2*data_for_mod[i*4+2]+data_for_mod[i*4+3]]);
		}
		break;
	case 64:
		for(int i=0;i<data_for_mod.size()/6;i++)
		{
			mod_out.push_back(_64QAM[32*data_for_mod[i*6]+16*data_for_mod[i*6+1]+8*data_for_mod[i*6+2]+4*data_for_mod[i*6+3]+2*data_for_mod[i*6+4]+data_for_mod[i*6+5]]);
		}
		break;
	default:
		break;
	}

}
void ProtocolE::repetition(int Burst_repetition,vector <complex <double> > &Mod_Out_temp,vector <complex <double> > &Rep_Out)
{
// 函数名: repetition
// 函数功能描述:对输入数据进行重复操作
// 输入参数: Burst_repetition、Mod_Out_temp
// 输出参数:Rep_Out
// Called By: run
// 补充说明:
// 修改日期: 2017/08/07：15:10:41
	vector <complex<double>>Rep_Out_temp;
	int k=Mod_Out_temp.size()/48;
	if(Burst_repetition==0)
	{
		Rep_Out.insert(Rep_Out.end(), Mod_Out_temp.begin(), Mod_Out_temp.end());
	}
		
	else
	{
		for(int i=0;i<k;i++)
		{
			for(int j=0;j<48;j++)
			Rep_Out_temp.push_back(Mod_Out_temp[j+i*48]);
			for (int k = 0; k <Burst_repetition; k++)
			{
				Rep_Out.insert(Rep_Out.end(), Rep_Out_temp.begin(), Rep_Out_temp.end());
			}
			Rep_Out_temp.clear();
		}
	}

	
}

void ProtocolE::PUSC()//括号内添加子信道选择模块参数！！
{
// 函数名: PUSC
// 函数功能描述:子载波排列规则，产生了子载波位置序列和导频位置序列
// 输入参数: 全局参数para.FFT_Size
// 输出参数:导频位置序列PiSubC，PiSubC_m和子载波位置序列DataSubC_Sym_New，DataSubC_Sym_New_m
// Called By: tx_param
// 补充说明:PUSC为部分使用的子信道，但使用全部的子信道也是可行的。本程序即使用了所有子信道
// 修改日期: 2017/08/07：16:13:45
    vector <int> UsedSubcarrier, PiCSet, PiVSet0, PiVSet1, PiSubC_m, DataSubC_Sym_T, DataSubC_Sym_mT, DataSubC_Sym,
    DataSubC_Sym_m, PermutationBase, Ps, subcarrier, SubCarrier_Seril;
	vector <int> PermutationBase12, PermutationBase8, PermutationBase6, PermutationBase5, PermutationBase4, PiCSet_m, DataSubC_Sym_LC, DataSubC_Sym_LC_m;
    int BitPerSymbolAftCode, DL_PermBase, N_SubChannel, N_SubCarrier, nk;
    switch (para.FFT_Size)
	{
	case 512:
		{
			int Renumbering_Seq_T[30] = {12, 13, 26, 9, 5, 15, 21, 6, 28, 4, 2, 7, 10, 18, 29, 17, 16, 3, 20, 24, 14, 8, 23, 1, 25, 27, 22, 19, 11, 0};
			vector <int> Renumbering_Seq(Renumbering_Seq_T,Renumbering_Seq_T+30);
			vector <int> PiCSet0, PiCSet1;
			//int PermutationBase[16] = { 6, 14, 2, 3, 10, 8, 11, 15, 9, 1, 13, 12, 5, 7, 4, 0 };
			PermutationBase5.push_back(4);
			PermutationBase5.push_back(2);
			PermutationBase5.push_back(3);
			PermutationBase5.push_back(1);
			PermutationBase5.push_back(0);

			for (int i = 47; i < 257; i++)
			{
				UsedSubcarrier.push_back(i);
			}
			for (int i = 258; i < 468; i++)
			{
				UsedSubcarrier.push_back(i);
			}

			PiCSet.push_back(4);
			PiCSet.push_back(8);

			PiCSet_m.push_back(0);
			PiCSet_m.push_back(12);


			for (int i = 1; i < 30; i++)
			{
				PiCSet.push_back(PiCSet[0]+14*i);
				PiCSet.push_back(PiCSet[1]+14*i);
				PiCSet_m.push_back(PiCSet_m[0]+14*i);
				PiCSet_m.push_back(PiCSet_m[1]+14*i);
			}


			para.PiSubC.insert(para.PiSubC.end(), PiCSet.begin(), PiCSet.end());

			for (int i = 0; i < para.PiSubC.size()/2; i++)
			{
				para.PiSubC[i] = para.PiSubC[i] + 47;
			}
			for (int i = para.PiSubC.size()/2; i < para.PiSubC.size(); i++)
			{
				para.PiSubC[i] = para.PiSubC[i] + 48;
			}
			sort(para.PiSubC.begin(), para.PiSubC.end());


			para.PiSubC_m.insert(para.PiSubC_m.end(), PiCSet_m.begin(), PiCSet_m.end());

			for (int i = 0; i < para.PiSubC_m.size()/2; i++)
			{
				para.PiSubC_m[i] = para.PiSubC_m[i] + 47;
			}
			for (int i = para.PiSubC_m.size()/2; i < para.PiSubC_m.size(); i++)
			{
				para.PiSubC_m[i] = para.PiSubC_m[i] + 48;
			}

			sort(para.PiSubC_m.begin(), para.PiSubC_m.end());
			para.LengthOfDataSubCarrier = UsedSubcarrier.size() - para.PiSubC.size();

			vector <int> Ismember_Flag;

			for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_T.push_back(i);
				for (int j = 0; j < para.PiSubC.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC[j])
					{
						Ismember_Flag.push_back(i);
						//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[76] - 76));
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
			{
				DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 46);
			}
			for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
			{
				DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 47);
			}
			for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_mT.push_back(i);
				for (int j = 0; j < para.PiSubC_m.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC_m[j])
					{
						Ismember_Flag.push_back(i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_mT.erase(DataSubC_Sym_mT.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
			{
				DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 46);
			}
			for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
			{
				DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 47);
			}
			DataSubC_Sym_LC.resize(360);
			DataSubC_Sym_LC_m.resize(360);
			for(int i = 0; i < Renumbering_Seq.size(); i++)  //！暂时只定义了first_DL_Zone的簇分配，从第二个Zone开始，簇分配需要参照协议再修改
			{
				copy(DataSubC_Sym.begin()+i*12,DataSubC_Sym.begin()+i*12+12,DataSubC_Sym_LC.begin()+12*Renumbering_Seq[i]);
				copy(DataSubC_Sym_m.begin()+i*12,DataSubC_Sym_m.begin()+i*12+12,DataSubC_Sym_LC_m.begin()+12*Renumbering_Seq[i]);
			}

			int N_Group = 3;
			int N_SubCarrier_Group = 120;

			N_SubChannel = 5; //子信道数
			N_SubCarrier = 24; //每个子信道数据子载波个数
			int subcarrier5[24][5];
			int DL_PermBase = 0;
			for (int k = 0; k < N_SubCarrier; k++)  //k
			{
				for (int s = 0; s < N_SubChannel; s++)  //s
				{
					nk = (k + 13 * s) % N_SubCarrier;
					circshift_ne(PermutationBase5, s, Ps);
					subcarrier5[k][s] = N_SubChannel*nk + (Ps[nk % N_SubChannel] + DL_PermBase) % N_SubChannel;
					//SubCarrier_Seril.push_back(subcarrier[k][s]);reshape的顺序不同。从列读入
					Ps.clear();
				}
			}
			for (int s = 0; s < N_SubChannel; s++)
			{
				for (int k = 0; k < N_SubCarrier; k++)
				{
					SubCarrier_Seril.push_back(subcarrier5[k][s]);
				}
			}

			for (int g = 1; g < N_Group; g++)
			{
				for (int i = 0; i < N_SubCarrier_Group; i++)
				{
					SubCarrier_Seril.push_back(SubCarrier_Seril[i]+N_SubCarrier_Group*g);
				}
			}


			//sort(SubCarrier_Seril.begin(), SubCarrier_Seril.end());
			for (int i = 0; i < SubCarrier_Seril.size(); i++)
			{
				para.DataSubC_Sym_New.push_back(DataSubC_Sym_LC[SubCarrier_Seril[i]]);
				para.DataSubC_Sym_m_New.push_back(DataSubC_Sym_LC_m[SubCarrier_Seril[i]]);
			}
			break;
		}
		
	case 1024:
		{
			int Renumbering_Seq_T[60] = {6, 48, 37, 21, 31, 40, 42, 56, 32, 47, 30, 33, 54, 18, 10, 15, 50, 51, 58, 46, 23, 45, 16, 57, 39, 35, 7, 55, 25, 59, 53, 11, 22, 38, 28, 19, 17, 3, 27, 12, 29, 26, 5, 41, 49, 44, 9, 8, 1, 13, 36, 14, 43, 2, 20, 24, 52, 4, 34, 0};
			vector <int> Renumbering_Seq(Renumbering_Seq_T,Renumbering_Seq_T+60);
		    vector <int> PiCSet0, PiCSet1;
            //int PermutationBase[16] = { 6, 14, 2, 3, 10, 8, 11, 15, 9, 1, 13, 12, 5, 7, 4, 0 };
            PermutationBase6.push_back(3);
			PermutationBase6.push_back(2);
			PermutationBase6.push_back(0);
			PermutationBase6.push_back(4);
			PermutationBase6.push_back(5);
			PermutationBase6.push_back(1);

			PermutationBase4.push_back(3);
			PermutationBase4.push_back(0);
			PermutationBase4.push_back(2);
			PermutationBase4.push_back(1);


		    for (int i = 93; i < 513; i++)
            {
                UsedSubcarrier.push_back(i);
            }
            for (int i = 514; i < 934; i++)
            {
                UsedSubcarrier.push_back(i);
            }

            PiCSet.push_back(4);
            PiCSet.push_back(8);

            PiCSet_m.push_back(0);
            PiCSet_m.push_back(12);


            for (int i = 1; i < 60; i++)
            {
                 PiCSet.push_back(PiCSet[0]+14*i);
				 PiCSet.push_back(PiCSet[1]+14*i);
				 PiCSet_m.push_back(PiCSet_m[0]+14*i);
				 PiCSet_m.push_back(PiCSet_m[1]+14*i);
            }

			
            para.PiSubC.insert(para.PiSubC.end(), PiCSet.begin(), PiCSet.end());

            for (int i = 0; i < para.PiSubC.size()/2; i++)
            {
                para.PiSubC[i] = para.PiSubC[i] + 93;
            }
			for (int i = para.PiSubC.size()/2; i < para.PiSubC.size(); i++)
            {
                para.PiSubC[i] = para.PiSubC[i] + 94;
            }
            sort(para.PiSubC.begin(), para.PiSubC.end());


            para.PiSubC_m.insert(para.PiSubC_m.end(), PiCSet_m.begin(), PiCSet_m.end());

            for (int i = 0; i < para.PiSubC_m.size()/2; i++)
            {
                para.PiSubC_m[i] = para.PiSubC_m[i] + 93;
            }
			for (int i = para.PiSubC_m.size()/2; i < para.PiSubC_m.size(); i++)
            {
                para.PiSubC_m[i] = para.PiSubC_m[i] + 94;
            }

			sort(para.PiSubC_m.begin(), para.PiSubC_m.end());
            para.LengthOfDataSubCarrier = UsedSubcarrier.size() - para.PiSubC.size();
           
			vector <int> Ismember_Flag;
	
            for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_T.push_back(i);
				for (int j = 0; j < para.PiSubC.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC[j])
					{
						Ismember_Flag.push_back(i);
						//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[76] - 76));
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
            {
                DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 92);
            }
            for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
            {
                DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 93);
            }
			for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_mT.push_back(i);
				for (int j = 0; j < para.PiSubC_m.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC_m[j])
					{
						Ismember_Flag.push_back(i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_mT.erase(DataSubC_Sym_mT.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
            {
                DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 92);
            }
            for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
            {
                DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 93);
            }

			DataSubC_Sym_LC.resize(720);
			DataSubC_Sym_LC_m.resize(720);
			for(int i = 0; i < Renumbering_Seq.size(); i++)  //！暂时只定义了first_DL_Zone的簇分配，从第二个Zone开始，簇分配需要参照协议再修改
			{
				//DataSubC_Sym_LC.replace(DataSubC_Sym_LC.begin()+12*Renumbering_Seq[i], 12,DataSubC_Sym,DataSubC_Sym.begin()+i, 12); 
				copy(DataSubC_Sym.begin()+i*12,DataSubC_Sym.begin()+i*12+12,DataSubC_Sym_LC.begin()+12*Renumbering_Seq[i]);
				copy(DataSubC_Sym_m.begin()+i*12,DataSubC_Sym_m.begin()+i*12+12,DataSubC_Sym_LC_m.begin()+12*Renumbering_Seq[i]);
				//DataSubC_Sym_LC_m.replace(DataSubC_Sym_LC_m.end(), DataSubC_Sym_m.begin()+(12*Renumbering_Seq[i]), DataSubC_Sym_m.begin()+(12*Renumbering_Seq[i])+12);
			}

			int N_Group = 3;
			int N_SubCarrier_Group = 240;

            N_SubChannel = 6; //子信道数
            N_SubCarrier = 24; //每个子信道数据子载波个数
            int subcarrier6[24][6];
			int DL_PermBase = 0;
            for (int k = 0; k < N_SubCarrier; k++)  //k
            {
                for (int s = 0; s < N_SubChannel; s++)  //s
                {
                    nk = (k + 13 * s) % N_SubCarrier;
                    circshift_ne(PermutationBase6, s, Ps);
                    subcarrier6[k][s] = N_SubChannel*nk + (Ps[nk % N_SubChannel] + DL_PermBase) % N_SubChannel ;
                    //SubCarrier_Seril.push_back(subcarrier[k][s]);reshape的顺序不同。从列读入
					Ps.clear();
                }
            }
			for (int s = 0; s < N_SubChannel; s++)
			{
				for (int k = 0; k < N_SubCarrier; k++)
				{
					SubCarrier_Seril.push_back(subcarrier6[k][s]);
				}
			}


			N_SubChannel = 4; //子信道数
            N_SubCarrier = 24; //每个子信道数据子载波个数
            int subcarrier4[24][4];
            for (int k = 0; k < N_SubCarrier; k++)  //k
            {
                for (int s = 0; s < N_SubChannel; s++)  //s
                {
                    nk = (k + 13 * s) % N_SubCarrier;
                    circshift_ne(PermutationBase4, s, Ps);
                    subcarrier4[k][s] = N_SubChannel*nk + (Ps[nk % N_SubChannel] + DL_PermBase) % N_SubChannel ;
                    //SubCarrier_Seril.push_back(subcarrier[k][s]);reshape的顺序不同。从列读入
					Ps.clear();
                }
            }
			for (int s = 0; s < N_SubChannel; s++)
			{
				for (int k = 0; k < N_SubCarrier; k++)
				{
					SubCarrier_Seril.push_back(144+subcarrier4[k][s]);
				}
			}

			for (int g = 1; g < N_Group; g++)
			{
				for (int i = 0; i < N_SubCarrier_Group; i++)
				{
					SubCarrier_Seril.push_back(SubCarrier_Seril[i]+N_SubCarrier_Group*g);
				}
			}


			//sort(SubCarrier_Seril.begin(), SubCarrier_Seril.end());
            for (int i = 0; i < SubCarrier_Seril.size(); i++)
            {
                para.DataSubC_Sym_New.push_back(DataSubC_Sym_LC[SubCarrier_Seril[i]]);
                para.DataSubC_Sym_m_New.push_back(DataSubC_Sym_LC_m[SubCarrier_Seril[i]]);
            }
            break;
        }
	
    case 2048:
		{
			int Renumbering_Seq_T[120] = {6, 108, 37, 81, 31, 100, 42, 116, 32, 107, 30, 93, 54, 78,
				10, 75, 50, 111, 58, 106, 23, 105, 16, 117, 39, 95, 7,
				115, 25, 119, 53, 71, 22, 98, 28, 79, 17, 63, 27, 72, 29,
				86, 5, 101, 49, 104, 9, 68, 1, 73, 36, 74, 43, 62, 20, 84,
				52, 64, 34, 60, 66, 48, 97, 21, 91, 40, 102, 56, 92, 47,
				90, 33, 114, 18, 70, 15, 110, 51, 118, 46, 83, 45, 76, 57,
				99, 35, 67, 55, 85, 59, 113, 11, 82, 38, 88, 19, 77, 3, 87,
				12, 89, 26, 65, 41, 109, 44, 69, 8, 61, 13, 96, 14, 103, 2,
				80, 24, 112, 4, 94, 0};
			vector <int> Renumbering_Seq(Renumbering_Seq_T,Renumbering_Seq_T+120);
			vector <int> PiCSet0, PiCSet1;
			//int PermutationBase[16] = { 6, 14, 2, 3, 10, 8, 11, 15, 9, 1, 13, 12, 5, 7, 4, 0 };
			PermutationBase12.push_back(6);
			PermutationBase12.push_back(9);
			PermutationBase12.push_back(4);
			PermutationBase12.push_back(8);
			PermutationBase12.push_back(10);
			PermutationBase12.push_back(11);
			PermutationBase12.push_back(5);
			PermutationBase12.push_back(2);
			PermutationBase12.push_back(7);
			PermutationBase12.push_back(3);
			PermutationBase12.push_back(1);
			PermutationBase12.push_back(0);

			PermutationBase8.push_back(7);
			PermutationBase8.push_back(4);
			PermutationBase8.push_back(0);
			PermutationBase8.push_back(2);
			PermutationBase8.push_back(1);
			PermutationBase8.push_back(5);
			PermutationBase8.push_back(3);
			PermutationBase8.push_back(6);


			for (int i = 185; i < 1025; i++)
			{
				UsedSubcarrier.push_back(i);
			}
			for (int i = 1026; i < 1866; i++)
			{
				UsedSubcarrier.push_back(i);
			}

			PiCSet.push_back(4);
			PiCSet.push_back(8);

			PiCSet_m.push_back(0);
			PiCSet_m.push_back(12);


			for (int i = 1; i < 120; i++)
			{
				PiCSet.push_back(PiCSet[0]+14*i);
				PiCSet.push_back(PiCSet[1]+14*i);
				PiCSet_m.push_back(PiCSet_m[0]+14*i);
				PiCSet_m.push_back(PiCSet_m[1]+14*i);
			}


			para.PiSubC.insert(para.PiSubC.end(), PiCSet.begin(), PiCSet.end());

			for (int i = 0; i < para.PiSubC.size()/2; i++)
			{
				para.PiSubC[i] = para.PiSubC[i] + 185;
			}
			for (int i = para.PiSubC.size()/2; i < para.PiSubC.size(); i++)
			{
				para.PiSubC[i] = para.PiSubC[i] + 186;
			}
			sort(para.PiSubC.begin(), para.PiSubC.end());


			para.PiSubC_m.insert(para.PiSubC_m.end(), PiCSet_m.begin(), PiCSet_m.end());

			for (int i = 0; i < para.PiSubC_m.size()/2; i++)
			{
				para.PiSubC_m[i] = para.PiSubC_m[i] + 185;
			}
			for (int i = para.PiSubC_m.size()/2; i < para.PiSubC_m.size(); i++)
			{
				para.PiSubC_m[i] = para.PiSubC_m[i] + 186;
			}

			sort(para.PiSubC_m.begin(), para.PiSubC_m.end());
			para.LengthOfDataSubCarrier = UsedSubcarrier.size() - para.PiSubC.size();

			vector <int> Ismember_Flag;

			for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_T.push_back(i);
				for (int j = 0; j < para.PiSubC.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC[j])
					{
						Ismember_Flag.push_back(i);
						//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			//DataSubC_Sym_T.erase(DataSubC_Sym_T.begin() + (Ismember_Flag[76] - 76));
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
			{
				DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 184);
			}
			for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
			{
				DataSubC_Sym.push_back(DataSubC_Sym_T[i] + 185);
			}
			for (int i = 0; i < UsedSubcarrier.size(); i++)
			{
				DataSubC_Sym_mT.push_back(i);
				for (int j = 0; j < para.PiSubC_m.size(); j++)
				{
					if (UsedSubcarrier[i] == para.PiSubC_m[j])
					{
						Ismember_Flag.push_back(i);
					}
				}
			}
			for (int i = 0; i < Ismember_Flag.size(); i++)
			{
				DataSubC_Sym_mT.erase(DataSubC_Sym_mT.begin() + (Ismember_Flag[i] - i));
			}
			Ismember_Flag.clear();
			for (int i = 0; i < para.LengthOfDataSubCarrier / 2; i++)
			{
				DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 184);
			}
			for (int i = para.LengthOfDataSubCarrier / 2; i < para.LengthOfDataSubCarrier; i++)
			{
				DataSubC_Sym_m.push_back(DataSubC_Sym_mT[i] + 185);
			}
			DataSubC_Sym_LC.resize(1440);
			DataSubC_Sym_LC_m.resize(1440);
			for(int i = 0; i < Renumbering_Seq.size(); i++)  //！暂时只定义了first_DL_Zone的簇分配，从第二个Zone开始，簇分配需要参照协议再修改
			{
				copy(DataSubC_Sym.begin()+i*12,DataSubC_Sym.begin()+i*12+12,DataSubC_Sym_LC.begin()+12*Renumbering_Seq[i]);
				copy(DataSubC_Sym_m.begin()+i*12,DataSubC_Sym_m.begin()+i*12+12,DataSubC_Sym_LC_m.begin()+12*Renumbering_Seq[i]);
			}

			int N_Group = 3;
			int N_SubCarrier_Group = 480;

			N_SubChannel = 12; //子信道数
			N_SubCarrier = 24; //每个子信道数据子载波个数
			int subcarrier12[24][12];
			int DL_PermBase = 0;
			for (int k = 0; k < N_SubCarrier; k++)  //k
			{
				for (int s = 0; s < N_SubChannel; s++)  //s
				{
					nk = (k + 13 * s) % N_SubCarrier;
					circshift_ne(PermutationBase12, s, Ps);
					subcarrier12[k][s] = N_SubChannel*nk + (Ps[nk % N_SubChannel] + DL_PermBase) % N_SubChannel + 1;
					//SubCarrier_Seril.push_back(subcarrier[k][s]);reshape的顺序不同。从列读入
					Ps.clear();
				}
			}
			for (int s = 0; s < N_SubChannel; s++)
			{
				for (int k = 0; k < N_SubCarrier; k++)
				{
					SubCarrier_Seril.push_back(subcarrier12[k][s]);
				}
			}


			N_SubChannel = 8; //子信道数
			N_SubCarrier = 24; //每个子信道数据子载波个数
			int subcarrier8[24][8];
			for (int k = 0; k < N_SubCarrier; k++)  //k
			{
				for (int s = 0; s < N_SubChannel; s++)  //s
				{
					nk = (k + 13 * s) % N_SubCarrier;
					circshift_ne(PermutationBase8, s, Ps);
					subcarrier8[k][s] = N_SubChannel*nk + (Ps[nk % N_SubChannel] + DL_PermBase) % N_SubChannel + 1;
					//SubCarrier_Seril.push_back(subcarrier[k][s]);reshape的顺序不同。从列读入
					Ps.clear();
				}
			}
			for (int s = 0; s < N_SubChannel; s++)
			{
				for (int k = 0; k < N_SubCarrier; k++)
				{
					SubCarrier_Seril.push_back(288+subcarrier8[k][s]);
				}
			}

			for (int g = 1; g < N_Group; g++)
			{
				for (int i = 0; i < N_SubCarrier_Group; i++)
				{
					SubCarrier_Seril.push_back(SubCarrier_Seril[i]+N_SubCarrier_Group*g);
				}
			}


			//sort(SubCarrier_Seril.begin(), SubCarrier_Seril.end());
			for (int i = 0; i < SubCarrier_Seril.size(); i++)
			{
				para.DataSubC_Sym_New.push_back(DataSubC_Sym_LC[SubCarrier_Seril[i] - 1]);
				para.DataSubC_Sym_m_New.push_back(DataSubC_Sym_LC_m[SubCarrier_Seril[i] - 1]);
			}
			break;
		}
}
}

void ProtocolE::IFFT_CP(vector<complex<double>> &IFFT_DataIn_All, vector<complex<double>> &Data_Tx)
{
// 函数名: IFFT_CP
// 函数功能描述:对将发送的调制数据进行IFFT操作以及加CP操作
// 输入参数: IFFT_DataIn_All
// 输出参数:Data_Tx
// Called By: run
// 补充说明:
// 修改日期: 2017/08/07：14:40:18
	Vector< complex<double> >  sn, xn;//注意，此处生成的是FFT库函数特定的Vector容器类，V是大写

	vector <complex <double> > IFFT_DataOut, IFFT_DataIn_Tmp, IFFT_DataIn_Tmp_1;

	for (int i = 0; i < IFFT_DataIn_All.size() / para.FFT_Size; i++)
	{
		IFFT_DataIn_Tmp.insert(IFFT_DataIn_Tmp.end(), IFFT_DataIn_All.begin() + i * para.FFT_Size, IFFT_DataIn_All.begin() + (i + 1) * para.FFT_Size);
		sn.resize(IFFT_DataIn_Tmp.size());
		for( int n=0; n<IFFT_DataIn_Tmp.size(); ++n )
			sn[n] = complex<double>( IFFT_DataIn_Tmp[n].real(), IFFT_DataIn_Tmp[n].imag() );//将已有数据转换到FFT库函数特定的Vector容器类对象sn中
		xn = ifftc2c( sn );//调用FFT库函数ifftc2c
		Vector< complex<double> >::const_iterator itrR = xn.begin();//定义FFT库函数特定的Vector容器类对象xn的迭代器itrR
		IFFT_DataIn_Tmp_1.resize(IFFT_DataIn_Tmp.size());
		for( int n=0; n<IFFT_DataIn_Tmp.size(); ++n )//将计算结果xn转换至输出容器IFFT_DataIn_Tmp中
		{
			IFFT_DataIn_Tmp[n].real((*itrR).real());
			IFFT_DataIn_Tmp[n].imag((*itrR++).imag());
		}

		for (int j = 0; j < para.FFT_Size; j++)
		{
			IFFT_DataIn_Tmp[j] = IFFT_DataIn_Tmp[j] *(double)sqrt((double)para.FFT_Size);//*sqrt(N)
		}
		IFFT_DataOut.insert(IFFT_DataOut.end(), IFFT_DataIn_Tmp.begin(), IFFT_DataIn_Tmp.end());
		IFFT_DataIn_Tmp.clear();
	}
	/**********************************************
						CP Insertion
	***********************************************/
	vector <complex <double> > CP_Out;
	int CP_Out_j = IFFT_DataOut.size() / para.FFT_Size;

	for (int i = 0; i  < CP_Out_j; i++)
	{
		for (int j = para.FFT_Size - para.LengthOfTgSp; j < para.FFT_Size; j++) //[896,1023]为重复部分，放在每列数据的头部
		{
			CP_Out.push_back(IFFT_DataOut[i * para.FFT_Size + j]);
		}
		for (int j = 0; j < para.FFT_Size; j++)
		{
			CP_Out.push_back(IFFT_DataOut[i * para.FFT_Size + j]);
		}

	}
	Data_Tx.insert(Data_Tx.end(), CP_Out.begin(), CP_Out.end());
}
				        
void ProtocolE::Clear()
{
	para.Preamble_out.clear();
	para.Poilts_Out.clear();

	para.DL_MAP.clear();
	para.MAC_head_DLMAP.clear();
	para.DL_MAP_out.clear();

	para.FCH_out.clear();
	para.DL_Frame_Prefix.clear();

	para.PiSubC.clear();
	para.PiSubC_m.clear();
	para.DataSubC_Sym_New.clear();
	para.DataSubC_Sym_m_New.clear();

	para.Next_SampleIndex.clear();
	para.Data_From_Pre.clear();
	para.I_Data.clear();
	para.Q_Data.clear();
	//para.Data_Tx.clear();
}
void ProtocolE::circshift(vector <int> &X, int k, vector <int> &Y)
{
// 函数名: circshift
// 函数功能描述:对输入数据进行循环右移k位
// 输入参数:X 
// 输出参数:Y
// Called By: randomization_generate、poilts_generate
// 补充说明:
// 修改日期: 2017/08/07：15:21:42
    vector <int> order;
    for (int i = X.size() - 1; i >= 0; i--)
    {
        order.push_back(i);
    }
	if (k == 0)//无移位操作
	{
		Y.insert(Y.end(), X.begin(), X.end());
	}
	else
	{
		for (int i = order[k] + 1; i < X.size(); i++)
		{
			Y.push_back(X[i]);
		}
		for (int i = 0; i <= order[k]; i++)
		{
			Y.push_back(X[i]);
		}
	}
}
void ProtocolE::circshift_ne(vector <int> &X, int k, vector <int> &Y)
{
// 函数名: circshift_ne
// 函数功能描述:对输入数据进行循环左移k位
// 输入参数: X
// 输出参数:Y
// Called By: PUSC
// 补充说明:
// 修改日期: 2017/08/07：15:23:26
    vector <int> order;
    for (int i = 0; i < X.size(); i++)
    {
        order.push_back(i);
    }
	if (k == 0)
	{
		Y.insert(Y.end(), X.begin(), X.end());
	}
	else
	{
		for (int i = order[k]; i < X.size(); i++)
		{
			Y.push_back(X[i]);
		}
		for (int i = 0; i < order[k]; i++)
		{
			Y.push_back(X[i]);
		}
	}
}
void ProtocolE::STBC_encode(vector <complex <double> > &Mod_In, int Nt, int NumOfUnit, vector <complex <double> > &Mod_Out,vector<complex<double>> &Mod_Out_2)
{
// 函数名: STBC_encode
// 函数功能描述:空时编码
// 输入参数: Mod_In、Nt、NumOfUnit、
// 输出参数:Mod_Out、Mod_Out_2
// Called By: run
// 补充说明:
// 修改日期: 2017/08/07：16:22:41
    int L = Mod_In.size() / NumOfUnit;//串并转换，对每一个子载波数据进行空时编码
    vector <complex<double>> temp,temp_2;
	/*ofstream mod_1,mod_2;
	mod_1.open("anten_1.txt",ofstream::out);
	mod_2.open("anten_2.txt",ofstream::out);*/
    switch (Nt)
    {
    case 1:
        Mod_Out.insert(Mod_Out.begin(), Mod_In.begin(), Mod_In.end());
        break;
    case 2:
        for(int i=0;i<NumOfUnit;i++)
		{
			for(int j=0;j<L;j=j+2)
			{
				temp.push_back(Mod_In[NumOfUnit*j+i]);//S1
				temp.push_back(-(conj(Mod_In[NumOfUnit*(j+1)+i])));//-conj(S2)
				temp_2.push_back(Mod_In[(j+1)*NumOfUnit+i]);//S2
				temp_2.push_back(conj(Mod_In[j*NumOfUnit+i]));//conj(S1)
			}
		}
		for(int i=0;i<L;i++)
		{
			for(int j=0;j<NumOfUnit;j++)
			{
				Mod_Out.push_back(temp[j*L+i]);
				Mod_Out_2.push_back(temp_2[j*L+i]);
				//mod_1<<temp[j*L+i]<<endl;
				//mod_2<<temp_2[j*L+i]<<endl;
			}
		}
		temp.clear();temp_2.clear();
        break;
    case 4:
        break;
    case 8:
        break;

    }
}

void ProtocolE::Subcarrier_Mapping(vector <complex <double> > &Mod_Out_Encoded,vector <complex <double> > &IFFT_DataIn)
{
// 函数名: Subcarrier_Mapping
// 函数功能描述:OFDMA数据子载波和导频子载波的数据映射
// 输入参数: Mod_Out_Encoded
// 输出参数:IFFT_DataIn
// Called By: run
// 补充说明:每个时隙分成两个24子载波。第一个24数据块映射到偶数OFDMA符号，第二个24数据块映射到奇数OFDMA符号。且每个子载波要乘以（1-2Wk）,具体值参考协议标准
// 修改日期: 2017/08/07：16:29:31
	vector <complex <double> > Mod_Out_Encoded_even, Mod_Out_Encoded_odd;
	for (int i=0;i<Mod_Out_Encoded.size()/48;i++)
	{
		for (int j=0;j<24;j++)
		{
			Mod_Out_Encoded_even.push_back(Mod_Out_Encoded[i*48+j]); //偶数符号的数据
		}
	}
	for (int i=0;i<Mod_Out_Encoded.size()/48;i++)
	{
		for (int j=0;j<24;j++)
		{
			Mod_Out_Encoded_odd.push_back(Mod_Out_Encoded[i*48+24+j]); //奇数符号的数据
		}
	}


	int poilts_off=92;//左侧保护带空闲子载波
	if(para.FFT_Size==512)
		poilts_off=46;
	else if(para.FFT_Size==2048)
		poilts_off=184;
	
	/*************************在所有符号中数据子载波插入数据*****************************/

	/******第一个，第二个符号数据*******/
	for(int j=0;j<para.NumOfUnit;j++)
	{
		IFFT_DataIn[para.DataSubC_Sym_New[j]]=Mod_Out_Encoded_even[j]*complex<double>(para.Poilts_Out[para.DataSubC_Sym_New[j]-poilts_off],0);//在使用的子信道中填充数据
	}

	for(int j=0;j<para.NumOfUnit;j++)
	{
		IFFT_DataIn[para.FFT_Size+para.DataSubC_Sym_m_New[j]]=Mod_Out_Encoded_odd[j]*complex<double>(para.Poilts_Out[para.DataSubC_Sym_m_New[j]-poilts_off+1],0);//在使用的子信道中填充数据
	}

	/*****第三个符号以后的数据********/

	/*****偶数符号*****/
	for(int i=1;i<(para.SymbolsForBurst2+2)/2;i++)
	{
		for(int j=para.Subchannel_offset_burst2*24;j<(para.Subchannel_offset_burst2+para.Subchannel_num_burst2)*24;j++)//从指定的子信道偏移处开始填充数据
		{
			IFFT_DataIn[2*i*para.FFT_Size+para.DataSubC_Sym_New[j]]=Mod_Out_Encoded_even[para.NumOfUnit+(i-1)*para.Subchannel_num_burst2*24+(j-para.Subchannel_offset_burst2*24)]*complex<double>(para.Poilts_Out[para.DataSubC_Sym_New[j]-poilts_off+2*i],0);//在使用的子信道中填充数据
		}               //2*i为符号偏移，j为逻辑子信道序列
	}
	/*****奇数符号*****/
	for(int i=1;i<(para.SymbolsForBurst2+2)/2;i++)
	{
		for(int j=para.Subchannel_offset_burst2*24;j<(para.Subchannel_offset_burst2+para.Subchannel_num_burst2)*24;j++)
		{
			IFFT_DataIn[(i*2+1)*para.FFT_Size+para.DataSubC_Sym_m_New[j]]=Mod_Out_Encoded_odd[para.NumOfUnit+(i-1)*para.Subchannel_num_burst2*24+(j-para.Subchannel_offset_burst2*24)]*complex<double>(para.Poilts_Out[para.DataSubC_Sym_m_New[j]-poilts_off+2*i+1],0);//在使用的子信道中填充数据
		}
	}

	/*************************在所有符号中插入导频序列*****************************/
	for (int i = 0; i < para.SymbolsForBurst/2; i++)
	{
		for (int j = 0; j < para.PiSubC.size(); j++)
		{
			IFFT_DataIn[para.PiSubC[j]-1+ 2 * i * para.FFT_Size] = para.Poilts_Out[para.PiSubC[j]-1-poilts_off+2*i];
		}
	}
	for (int i = 0; i < para.SymbolsForBurst/2; i++)
	{
		for (int j = 0; j < para.PiSubC_m.size(); j++)
		{
			IFFT_DataIn[para.PiSubC_m[j]-1+ (2 * i+1) * para.FFT_Size] = para.Poilts_Out[para.PiSubC_m[j]-1-poilts_off+2*i+1];
		}
	}

}

int ProtocolE::rand_max()
{
	double rand_num;
	srand((unsigned)time(0));
	return rand_num = rand();
}
double ProtocolE::randn()
{
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1)
        return randn();
    double c = sqrt(-2 * log(r) / r);
	double y = u * c;
	/*if (y >= 0)
	{
		y = y - 1;
	}*/
	if (y <= 0)
	{
		y = -y;
	}
	y = (int)y % 2;
    return y;
}
void ProtocolE::de2bi(int n_decimal, int bin_size, vector <int> &a_binary)
{
    int i_binary = 0;
    if (n_decimal <= 0)
    {
        for (int i = 0; i < bin_size; i++)
        {
            a_binary.push_back(0);
        }
    }
    else
    {
        a_binary.push_back(n_decimal % 2);
        while(n_decimal != 0 && n_decimal != 1)
        {
            n_decimal = n_decimal / 2;
            i_binary = n_decimal % 2;
            a_binary.push_back(i_binary);
        }
		while (a_binary.size() != bin_size)
		{
			a_binary.push_back(0);
		}
    }
}
int ProtocolE::c2i(char ch)
{
    // 如果是数字，则用数字的ASCII码减去48, 如果ch = '2' ,则 '2' - 48 = 2
    if(isdigit(ch))
        return ch - 48;

    // 如果是字母，但不是A~F,a~f则返回
    if( ch < 'A' || (ch > 'F' && ch < 'a') || ch > 'z' )
        return -1;

    // 如果是大写字母，则用数字的ASCII码减去55, 如果ch = 'A' ,则 'A' - 55 = 10
    // 如果是小写字母，则用数字的ASCII码减去87, 如果ch = 'a' ,则 'a' - 87 = 10
    if(isalpha(ch))
        return isupper(ch) ? ch - 55 : ch - 87;

    return -1;
}
void ProtocolE::hex2bit(char *hex, vector <int> &b_binary)
{
/*功能：将十六进制字符串转换为整型(int)数值*/
	int len;
	vector <int> tmp_binary;
    int temp = 0;

    // 此例中 hex = "1de" 长度为3, hex是main函数传递的
    len = strlen(hex);

    for (int i = 0; i < len; i++)
	{
		temp = c2i(*(hex + i));
		switch (temp)
		{
		case 0:
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			break;
		case 1:
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			break;
		case 2:
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			break;
		case 3:
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			break;
		case 4:
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			break;
		case 5:
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			break;
		case 6:
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			break;
		case 7:
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			break;
		case 8:
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			break;
		case 9:
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			break;
		case 10:
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			break;
		case 11:
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			break;
		case 12:
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(0);
			break;
		case 13:
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			tmp_binary.push_back(1);
			break;
		case 14:
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(0);
			break;
		case 15:
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			tmp_binary.push_back(1);
			break;
		}
		b_binary.insert(b_binary.end(), tmp_binary.begin(), tmp_binary.end());
		tmp_binary.clear();
	}
}
void ProtocolE::bitrp (vector <complex <double> > &Data_In, int n)
{
	// 位反转置换 Bit-reversal Permutation
	int i, j, a, b, p;
	complex<double> temp;
	for (i = 1, p = 0; i < n; i *= 2)
	{
		p ++;
	}
	for (i = 0; i < n; i ++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j ++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if ( b > i)
		{
			temp=Data_In[b];
			Data_In[b]=Data_In[i];
			Data_In[i]=temp;
		}
	}
}
void ProtocolE::de2bi_1(int deData[], int len, int biData[])
{
	for (int j = 0; j<len; j++)
	{
		for (int i = 0; i <= 7; i++)
		{
			int t = 1 << 7 - i;
			int k = deData[j] & t;
			if (k)
				biData[j * 8 + i] = 1;
			else
				biData[j * 8 + i] = 0;
		}
	}
}
void ProtocolE::bi2de_1(int biData[], int len, int deData[])
{
	double a = 2;
	for (int j = 0; j<len / 8; j++)
	{
		for (int i = 0; i <8; i++)	
			deData[j]+= biData[j * 8 + i] * pow(a, 7 - i);
	}
}
void ProtocolE::subblock_interleaver(int*Data_In, int len,int M, int J, int* Data_Out)
{ //subblock_interleaver
	int ii, jj, tmp, T;
	//char T;
	int k = 0;
	int p = (int)pow((double)2, M);
	int f, m;
	int* BRO = new int[p];
	int* BRO_1 = new int[p];
	int* biBRO = new int[8 * p];
	for (ii = 0; ii<p; ii++) { //for
		BRO[ii] = ii;
		BRO_1[ii]=0;
	} 
	de2bi_1(BRO, p, biBRO);

	for (ii = 0; ii<8*p; ii += 8) { //for

		for (jj = 0; jj<M/2; jj++) { //for
			tmp=biBRO[ii+jj+8-M];
			biBRO[ii + jj + 8 - M] = biBRO[ii + 7 - jj];
			biBRO[ii + 7 - jj] = tmp;
		}
	} //for

	bi2de_1(biBRO,8*p ,BRO_1);

	for (ii = 0; ii < len / 2; ii++)
	{
		f = (int)floor((float)(k / J));
		m = (int)fmod((float)k, (float)J);
		T = p*(m)+BRO_1[f];
		//T = 2 ^ M*(mod(k, J)) + bitrevorder(floor(k / J) + 1);
		while (T >= len / 2)
		{
			k = k + 1;
			f = (int)floor(float(k / J));
			m = (int)fmod((float)k, (float)J);
			T = p*(m) + BRO_1[f];
			//T = 2 ^ M*(mod(k, J)) + bitrevorder(floor(k / J) + 1);
		}

		Data_Out[ii * 2] = Data_In[T * 2];
		Data_Out[ii * 2 + 1] = Data_In[T * 2 + 1];
		k++;
	}
}

void ProtocolE::ctc_interleaver(int*Data_In, int len,
	char P0, char P1, char P2, char P3, int* Data_Out)
{
	short  ii, jj;
	char   tmp1;
	char   data[864]; //QAM64:27*16*2 
	//short  P1 = 3 * datain.size / 8;
	//Step 1: Swich alternate Couples
	for (jj = 0; jj<len; jj += 4) { //for
		tmp1 = Data_In[jj + 2];
		Data_In[jj + 2] = Data_In[jj + 3];
		Data_In[jj + 3] = tmp1;
	}//for

	//Step 2: Pi(j)
	for (jj = 0; jj < (len / 2); jj++) { //for
		ii = ((jj % 4) == 1) ? (P0*jj + 1 + (len / 4) + P1) % (len / 2) :
			((jj % 4) == 2) ? (P0*jj + 1 + P2) % (len / 2) :
			((jj % 4) == 3) ? (P0*jj + 1 + (len / 4) + P3) % (len / 2) :
			(P0*jj + 1) % (len / 2);

	data[jj * 2 + 0] = Data_In[ii * 2 + 0];
	data[jj * 2 + 1] = Data_In[ii * 2 + 1];
	//printf("interleaver (%d=>%d) 0b%x%x\n", jj, ii, data[ii*2-2] , data[ii*2-1]  ) ;
	} //for
	//copy back to datain
	for (int i = 0; i<len; i++) { //for
		Data_Out[i] = data[i];
	} //for
} //ctc_interleaver
void  ProtocolE::constituent_encoder(char aa, char bb, char *state, int *out1, int *out2) {
	char  tmp0, tmp1, tmp2;

	//tmp0 = (*state & 0x1) ^ ((*state >> 2) & 0x1) ^ aa ^ bb;
	//tmp0 += (((*state & 0x1) ^ bb) & 0x1) << 1;
	//tmp0 += ((((*state >> 1) & 0x1) ^ bb) & 0x1) << 2;
	//tmp1 = ((*state >> 2) & 0x1) ^ ((*state >> 1) & 0x1) ^ (tmp0 & 0x1);
	//tmp2 = ((*state >> 2) & 0x1) ^ (tmp0 & 0x1);
	//*state = tmp0 & 0x7;
	//*out1 = tmp1;
	//*out2 = tmp2;

	tmp0 = ((((*state >> 1) & 0x1) ^ bb) & 0x1);
	tmp0 += ((((*state >> 2) & 0x1) ^ bb) & 0x1) << 1;
	tmp0 += ((*state & 0x1) ^ ((*state >> 2) & 0x1) ^ aa ^ bb) << 2;

	tmp1 = (*state & 0x1) ^ ((*state >> 1) & 0x1) ^ ((tmp0>>2) & 0x1);
	tmp2 = (*state & 0x1) ^ ((tmp0>>2) & 0x1);
	*state = tmp0 & 0x7;
	*out1 = tmp1;
	*out2 = tmp2;
} //constituent_encoder
float ProtocolE::round(float data)
{
//	int data_int = (float) floor((float) data);
//	int data_fo = data - floor(data);
	if (data - floor(data) > 0.5)
	{
		data = floor(data) + 1;
	}
	else
	{
		data = floor(data);
	}
	return data;
}
float ProtocolE::rand1()
{
	int i;
	//srand((unsigned)time(0));
	double temp;
	temp = (double)rand()/( RAND_MAX + 1 );
	//temp = ( int )( temp * 2 );
	if (temp > 0.9)
	{
		temp = temp - 0.5;
	}
	return temp;
}




void ProtocolE::channle_init_996(float LengthOfBurst, float IntervalBetweenFrames, float MobileSpeed, float NumOfFreq, int NumOfChannels, float Channel_Type, float uplink,
								vector <float> &Path_Delay, vector <float>&Path_Average_Atmp, vector <float>&aglT, vector <float>&An, vector <float>&Bn, vector <float>&Wn, vector <float>&RndPhase, vector <float>&aglR,
								vector <int> &distypeR, vector <int> &distypeT,
								int *UpdateInterval, int *UpdatesPerBurst, float *dt, int *NumOfTaps, 
								int *aglsprR, int *aglsprT, float *dR, float *dT, 
								int *corrmodel, int *alf, float *fd, int *loson)
{
	//vector <float> Path_Delay, Path_Average_Atmp, aglT, An, Bn, Wn, RndPhase, aglR;
	//vector <int> distypeR, distypeT;
	/*vector <float> Path_Average_Atmp;
	vector <float> aglT;
	vector <float> An;
	vector <float> Bn;
	vector <float> Wn;
	vector <float> RndPhase;
	vector <float> aglR;*/
    *dt = IntervalBetweenFrames / LengthOfBurst;
    int MaxPhaseChange, M_OSCIL, M_TONES, N_MAGIC;
    float Path_Power_sum, fd_max, LastUpdateInterval;
    double CarrierFrequency, SpeedOfLight;
    vector <float> Path_Power_dB, Path_Power, belta;
    if (Channel_Type > 1.09 && Channel_Type < 1.11)//判断Channel_Type==1.1
    {
        *NumOfTaps = 5;

        Path_Delay.push_back(round(0 * 1e-9 / *dt));
        Path_Delay.push_back(round(0 * 1e-9 / *dt));
        Path_Delay.push_back(round(110 * 1e-9 / *dt));
        Path_Delay.push_back(round(190 * 1e-9 / *dt));
        Path_Delay.push_back(round(410 * 1e-9 / *dt));

        Path_Power_dB.push_back(0);
        Path_Power_dB.push_back(-6.51);
        Path_Power_dB.push_back(-16.21);
        Path_Power_dB.push_back(-25.71);
        Path_Power_dB.push_back(-29.31);

        for (int i = 0; i < Path_Power_dB.size(); i++)
        {
            Path_Power.push_back(pow(10, Path_Power_dB[i] / 10));
        }

        Path_Power_sum = Path_Power[1];
        for (int i = 1; i < Path_Power.size(); i++)
        {
            Path_Power_sum = Path_Power_sum + Path_Power[i];
        }
        for (int i = 0; i < Path_Power.size(); i++)
        {
            Path_Average_Atmp.push_back( sqrt(Path_Power[i] / Path_Power_sum) );
        }
        *corrmodel = 1;
        *alf = 0;
        *loson = 1;

        if (uplink == 1)
        {
            aglT.push_back(22.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);

            *aglsprT = 35;

            distypeT.push_back(1);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT = 0.5;

            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);

            *aglsprR = 2;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 4;
        }
        else
        {
            aglR.push_back(22.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);

            *aglsprR = 35;

            distypeR.push_back(1);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 0.5;

            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);

            *aglsprT = 2;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT = 4;
        }
    }

    else if (Channel_Type > 1.19 && Channel_Type < 1.21)//判断Channel_Type==1.2
    {
        *NumOfTaps = 4;

        Path_Delay.push_back( round( 0 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 110 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 190 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 410 * 1e-9 / *dt) );

        Path_Power_dB.push_back(0);
        Path_Power_dB.push_back(-9.7);
        Path_Power_dB.push_back(-19.2);
        Path_Power_dB.push_back(-22.8);

        for (int i = 0; i < Path_Power_dB.size(); i++)
        {
            Path_Power.push_back( pow(10, Path_Power_dB[i] / 10) );
        }

        Path_Power_sum = Path_Power[0];
        for (int i = 1; i < Path_Power.size(); i++)
        {
            Path_Power_sum = Path_Power_sum + Path_Power[i];
        }
        for (int i = 0; i < Path_Power.size(); i++)
        {
            Path_Average_Atmp.push_back( Path_Power[i] / Path_Power_sum);
        }

        *corrmodel = 1;
        *alf = 0;
        *loson = 0;

        if (uplink ==1)
        {
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);

            *aglsprT = 35;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT =0.5;

            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);

            *aglsprR = 2;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 4;
        }
        else
        {
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);

            *aglsprR = 35;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 0.5;

            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);

            *aglsprT = 2;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT = 4;
        }
    }

    else if (Channel_Type == 2)
    {
        *NumOfTaps = 6;

        Path_Delay.push_back( round( 0 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 310 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 710 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 1090 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 1730 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 2510 * 1e-9 / *dt) );

        Path_Power_dB.push_back(0);
        Path_Power_dB.push_back(-1);
        Path_Power_dB.push_back(-9);
        Path_Power_dB.push_back(-10);
        Path_Power_dB.push_back(-15);
        Path_Power_dB.push_back(-20);

        for (int i = 0; i < Path_Power_dB.size(); i++)
        {
            Path_Power.push_back( pow(10, Path_Power_dB[i] / 10) );
        }

        Path_Power_sum = Path_Power[0];
        for (int i = 1; i < Path_Power.size(); i++)
        {
            Path_Power_sum = Path_Power_sum + Path_Power[i];
        }
        for (int i = 0; i < Path_Power.size(); i++)
        {
            Path_Average_Atmp.push_back( Path_Power[i] / Path_Power_sum);
        }

        *corrmodel = 1;
        *alf = 22.5;
        *loson = 0;

        if (uplink == 1)
        {
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);
            aglT.push_back(67.5);

            *aglsprT = 35;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT =0.5;

            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);

            *aglsprR = 2;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 4;
        }
        else
        {
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);
            aglR.push_back(67.5);

            *aglsprR = 35;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR =0.5;

            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);

            *aglsprT = 2;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT = 4;
        }
    }

    else if (Channel_Type == 3)
    {
        *NumOfTaps = 6;

        Path_Delay.push_back( round( 0 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 200 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 800 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 1200 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 2300 * 1e-9 / *dt) );
        Path_Delay.push_back( round( 3700 * 1e-9 / *dt) );

        Path_Power_dB.push_back(0);
        Path_Power_dB.push_back(-0.9);
        Path_Power_dB.push_back(-4.9);
        Path_Power_dB.push_back(-8);
        Path_Power_dB.push_back(-7.8);
        Path_Power_dB.push_back(-23.9);

        for (int i = 0; i < Path_Power_dB.size(); i++)
        {
            Path_Power.push_back( pow(10, Path_Power_dB[i] / 10) );
        }

        Path_Power_sum = Path_Power[0];
        for (int i = 1; i < Path_Power.size(); i++)
        {
            Path_Power_sum = Path_Power_sum + Path_Power[i];
        }
        for (int i = 0; i < Path_Power.size(); i++)
        {
            Path_Average_Atmp.push_back( Path_Power[i] / Path_Power_sum);
        }

        *corrmodel = 1;
        *alf = 22.5;
        *loson = 0;

        if (uplink == 1)
        {
            aglT.push_back(22.5);
            aglT.push_back(67.5);
            aglT.push_back(22.5);
            aglT.push_back(67.5);
            aglT.push_back(22.5);
            aglT.push_back(67.5);

            *aglsprT = 35;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT =0.5;

            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);
            aglR.push_back(50);

            *aglsprR = 2;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR = 4;
        }
        else
        {
            aglR.push_back(22.5);
            aglR.push_back(67.5);
            aglR.push_back(22.5);
            aglR.push_back(67.5);
            aglR.push_back(22.5);
            aglR.push_back(67.5);

            *aglsprR = 35;

            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);
            distypeR.push_back(2);

            *dR =0.5;

            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);
            aglT.push_back(50);

            *aglsprT = 2;

            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);
            distypeT.push_back(2);

            *dT = 4;
        }
    }

    else
    {
        *NumOfTaps = 1;
        Path_Delay.push_back(0);
        Path_Power_dB.push_back(0);
        Path_Power.push_back(1);
        Path_Average_Atmp.push_back(1);

        *corrmodel = 0;
        *alf = 0;
        *loson = 0;

        aglT.push_back(0);
        *aglsprT = 104;
        distypeT.push_back(2);
        *dT =0.5;
        aglR.push_back(50);
        *aglsprR = 2;
        distypeR.push_back(2);
        *dR = 4;
    }

    CarrierFrequency = 2.3e9;
    SpeedOfLight = 3e8;
    fd_max = MobileSpeed / 3.6 / SpeedOfLight * CarrierFrequency;
    MaxPhaseChange = 5;

    if (MobileSpeed != 0)
    {
        fd_max = MobileSpeed / 3.6 / SpeedOfLight * CarrierFrequency * cos(*alf * 2 * pi / 360);

        *UpdateInterval = round((MaxPhaseChange / 360) / fd_max / *dt / 1000);

//****************限定"UpdateInterval" 在区间[1;LengthOfburst]****************
        if (*UpdateInterval > 1)
        {
            *UpdateInterval = *UpdateInterval;
        }
        else if (*UpdateInterval < 1)
        {
            *UpdateInterval = 1;
        }
        if (*UpdateInterval > LengthOfBurst)
        {
            *UpdateInterval = LengthOfBurst;
        }
        else if (*UpdateInterval < LengthOfBurst)
        {
            *UpdateInterval = *UpdateInterval;
        }

        *UpdatesPerBurst = round((LengthOfBurst - 1) / *UpdateInterval); // 每帧中信道参数更新次数
        LastUpdateInterval = LengthOfBurst - *UpdatesPerBurst * *UpdateInterval; // 更新间隔
    }
    else
    {
        *UpdatesPerBurst = 0;
        *UpdateInterval = 0;
        LastUpdateInterval = LengthOfBurst;
    }

//******************* Parameters for fading generating ***************************
    M_OSCIL = NumOfFreq;
    M_TONES = M_OSCIL + 1;
    N_MAGIC = M_OSCIL * 4 + 2;

//******************* Initialise the coefficients for fading generating *****************
    for (int i = 0; i < M_OSCIL; i++)
    {
        belta.push_back(pi * i / M_OSCIL);
        An.push_back(2 * cos(belta[i]));
        Bn.push_back(2 * sin(belta[i]));
    }
    An.push_back(1);
    Bn.push_back(1);

//******************* Initialise the tone frequencies *****************
    for (int i = 0; i < M_OSCIL; i++)
    {
        Wn.push_back(2 * pi * fd_max * cos(2 * pi * i / N_MAGIC));
    }
    Wn.push_back(2 * pi * fd_max);

	RndPhase.clear();
	srand((unsigned)time(0));
    for (int i = 0; i < NumOfChannels; i++)
    {
        for (int j = 0; j < *NumOfTaps; j++)
        {
            for(int k = 0; k < M_TONES; k++)
			{
				RndPhase.push_back(1 - 2 * rand1());
			}

        }
    }
	*fd = fd_max;
}

void ProtocolE::JakesIV(double Path_Ampli,double *An,double *Bn,double *Wn,double *Phase,int NumOfFreq,int UpdateInterval,\
						int UpdatesPerBurst,int LengthOfBurst,int t_tmp,double dt,double *fading_re,double *fading_im)
{  
   int i,j,k,n,M_tones,M_magic;
   double cosine, scale,ival,qval;
   

   M_tones=NumOfFreq+1;
   M_magic=NumOfFreq*4+2;//N//

   i=0;
   scale=1.414*Path_Ampli/sqrt((float)M_magic);
   
   for(k=0;k<=UpdatesPerBurst;k++)
   {
      //Reset sum 
      ival = 0;
      qval = 0;
      for(n=0;n<M_tones;n++)
      {
	cosine = cos(Wn[n]*t_tmp*dt+Phase[n]);
		
        ival = ival+An[n] * cosine;
        qval = qval+Bn[n] * cosine;
      }
      
      for(j=i;j<i+UpdateInterval;j++)
      {
      	  /*if(j<=UpdateInterval*UpdatesPerBurst)*/
      	  if(j<LengthOfBurst)
      	  {
             *(fading_re+j)= ival * scale;
             *(fading_im+j)= qval * scale;
          }
        
      }
      i=i+UpdateInterval;
      t_tmp=t_tmp+UpdateInterval;    
   } 
}

void ProtocolE::mimocorr_cx(double anglespread,double angle,double d,int M,int distype,double *correlation2_r,double *correlation2_i,int BS)
{
	int L=1000;
    double anglespread1=180;
	double anglerange=180;//
    double c=0;//
	int ii,m,n,jj;

	double *p;//
	double *fai1;//
	double *FAI;//
	double *matrix_r;
	double *matrix_i;

	p=(double *)calloc(8,L);
	fai1=(double *)calloc(8,L);
	FAI=(double *)calloc(8,L);
	matrix_r=(double *)calloc(8,M);
	matrix_i=(double *)calloc(8,M);

	if(anglespread==0)
	{
		anglespread=0.00000001;
	}

	for(ii=0;ii<L;ii++)
	{
		if(distype==1)//uniform distribution
		{
			fai1[ii]=angle-anglerange+2*anglerange*(ii+1)/L;
     		FAI[ii]=d*sin(2*pi*fai1[ii]/360);
    		p[ii]=1.0/L;
	    	//c=c+p[ii];
		}
		else if(distype==2)//Laplace distribution
		{
			fai1[ii]=angle-anglespread1+2*anglespread1*(ii+1)/L;
			FAI[ii]=d*sin(2*pi*fai1[ii]/360);
			p[ii]=sqrt((float)2)*anglespread1*exp(-sqrt((float)2)*fabs(fai1[ii]-angle)/anglespread)/(anglespread*L);//Laplace?
			//c=c+p[ii];
		}
		else//normal distribution
		{
			//mexErrMsgTxt("distribution type error");
			fai1[ii]=angle-anglespread1+2*anglespread1*(ii+1)/L;
			FAI[ii]=d*sin(2*pi*fai1[ii]/360);
			p[ii]=exp(-(fai1[ii]-angle)*(fai1[ii]-angle)/(2*anglespread*anglespread))/(anglespread*L);
			//c=c+p[ii];
		}
		if (BS==1)
		{
			if(fai1[ii]>180)
			{
	    		fai1[ii]=fai1[ii]-360;
			}
	    	if(fai1[ii]<-180)
			{
    			fai1[ii]=fai1[ii]+360;
			}
    		if(fabs(fai1[ii])>90)
			{
	    		p[ii]=p[ii]/pow((float)10,(float)20/20);
			}
	    	else
			{
	    		p[ii]=p[ii]/pow(10,fai1[ii]*fai1[ii]/4900*0.6);
			}
		}
		c=c+p[ii];
	}

	for(ii=0;ii<M*M;ii++)
	{
		correlation2_r[ii]=0;
		correlation2_i[ii]=0;
	}

	for(m=0;m<L;m++)
	{
		for(n=0;n<M;n++)
		{
			matrix_r[n]=cos(FAI[m]*2*pi*n);
			matrix_i[n]=-sin(FAI[m]*2*pi*n);
		}
		for(ii=0;ii<M;ii++)
		{
			for(jj=0;jj<M;jj++)
			{
				correlation2_r[ii*M+jj]+=(matrix_r[ii]*matrix_r[jj] + matrix_i[ii]*matrix_i[jj]) * p[m]/c;
				correlation2_i[ii*M+jj]+=( - matrix_r[ii]*matrix_i[jj] + matrix_i[ii]*matrix_r[jj]) * p[m]/c;
			}
		}
	}
	free(p);
	free(fai1);
	free(FAI);
	free(matrix_r);
	free(matrix_i);
}

int ProtocolE::cx_chol(double* a_r,double* a_i,int n)
{
	int i,j,k,u,l;
    if ((a_r[0]+1.0==1.0)||(a_r[0]<0.0))
    {
		return(-2);
	}
    a_r[0]=sqrt(a_r[0]);
    for (i=1; i<=n-1; i++)
    {
		u=i*n; 
		a_r[u]=a_r[u]/a_r[0];
		a_i[u]=a_i[u]/a_r[0];
	}
    for (j=1; j<=n-1; j++)
    {
		l=j*n+j;
        for (k=0; k<=j-1; k++)
        { 
			u=j*n+k; 
			//a_r[l]=a_r[l]-a_r[u]*a_r[u];
			a_r[l]=a_r[l]-a_r[u]*a_r[u]-a_i[u]*a_i[u];//(1)
		}
        if ((a_r[l]+1.0==1.0)||(a_r[l]<0.0))
        {
			return(-2);
		}
        a_r[l]=sqrt(a_r[l]);
        //d=d*a[l];
        for (i=j+1; i<=n-1; i++)
        {
			u=i*n+j;
            for (k=0; k<=j-1; k++)
			{
				//a_r[u]=a_r[u]-a_r[i*n+k]*a_r[j*n+k];
				a_r[u]=a_r[u]-a_r[i*n+k]*a_r[j*n+k]-a_i[i*n+k]*a_i[j*n+k];
				a_i[u]=a_i[u]+a_r[i*n+k]*a_i[j*n+k]-a_i[i*n+k]*a_r[j*n+k];
			}
            a_r[u]=a_r[u]/a_r[l];
			a_i[u]=a_i[u]/a_r[l];//(2)
        }
    }
    for (i=0; i<=n-2; i++)
	{
		for (j=i+1; j<=n-1; j++)
		{
           a_r[i*n+j]=0.0;
		   a_i[i*n+j]=0.0;
		}
	}
    return(2);
}

void ProtocolE::mimochannel_cx(int Nr,int Nt,double aglsprR,double aglR,double aglsprT,double aglT,double dR,double dT,int distypeR,int distypeT,double *hr_r,double *hr_i,int uplink)
{
	double *corrR_r;
	double *corrR_i;
	double *corrT_r;
	double *corrT_i;
	int ri,rj;
	int ti,tj;

	corrR_r = (double*)calloc(8,Nr*Nr);//强行更变void到double，试可行否？
	corrR_i = (double*)calloc(8,Nr*Nr);
	corrT_r = (double*)calloc(8,Nt*Nt);
	corrT_i = (double*)calloc(8,Nt*Nt);
	//mexErrMsgTxt("mimochannel cx");
	if (uplink==1)
	{
		mimocorr_cx(aglsprR,aglR,dR,Nr,distypeR,corrR_r,corrR_i,1);
		mimocorr_cx(aglsprT,aglT,dT,Nt,distypeT,corrT_r,corrT_i,0);
	}
	else //downlink
	{
		mimocorr_cx(aglsprR,aglR,dR,Nr,distypeR,corrR_r,corrR_i,0);
		mimocorr_cx(aglsprT,aglT,dT,Nt,distypeT,corrT_r,corrT_i,1);
	}

    //kron
	for(ri=0;ri<Nr;ri++)
	{
		for(rj=0;rj<Nr;rj++)
		{
			for(ti=0;ti<Nt;ti++)
			{
				for(tj=0;tj<Nt;tj++)
				{
					//hr_r[(ri*Nr+ti)*Nr*Nt+(rj*Nr+tj)]=corrR_r[ri*Nr+rj]*corrT_r[ti*Nt+tj]-corrR_i[ri*Nr+rj]*corrT_i[ti*Nt+tj];
					//hr_i[(ri*Nr+ti)*Nr*Nt+(rj*Nr+tj)]=corrR_r[ri*Nr+rj]*corrT_i[ti*Nt+tj]+corrR_i[ri*Nr+rj]*corrT_r[ti*Nt+tj];
					hr_r[(ri*Nt+ti)*Nr*Nt+(rj*Nt+tj)]=corrR_r[ri*Nr+rj]*corrT_r[ti*Nt+tj]-corrR_i[ri*Nr+rj]*corrT_i[ti*Nt+tj];
					hr_i[(ri*Nt+ti)*Nr*Nt+(rj*Nt+tj)]=corrR_r[ri*Nr+rj]*corrT_i[ti*Nt+tj]+corrR_i[ri*Nr+rj]*corrT_r[ti*Nt+tj];
				}
			}
		}
	}
	//bchol(hr,Nr*Nt,&det);
	cx_chol(hr_r,hr_i,Nr*Nt);

	free(corrR_r);
	free(corrR_i);
	free(corrT_r);
	free(corrT_i);
}

void ProtocolE::mul_cx(double *hr_r,double *hr_i,double *b_r,double *b_i,int m,int n,int k,double *c_r,double *c_i)
//hr stored in row,b stored in column,c stored in column
{
	int i,j,l,u;
    for (i=0; i<=m-1; i++)
	{
		for (j=0; j<=k-1; j++)
		{
			u=i+j*m; 
			c_r[u]=0.0;
			c_i[u]=0.0;
			for (l=0; l<=n-1; l++)
			{
				c_r[u]=c_r[u]+hr_r[i*n+l]*b_r[l+j*n]-hr_i[i*n+l]*b_i[l+j*n];
				c_i[u]=c_i[u]+hr_r[i*n+l]*b_i[l+j*n]+hr_i[i*n+l]*b_r[l+j*n];
			}
		}
	}
}

void ProtocolE::corrj_channel(	double *InData_pr, double *InData_pi, double *An, double *Bn, double *Wn, \
			double *RndPhase, int NumOfFreq, int UpdateInterval, int UpdatesPerBurst, \
			long LengthOfBurst, int SampleIndex, double dt, int NumOfTaps, \
			double *Path_Delay, double *Path_Average_Amp,  \
			double *fore_data_pr,double *fore_data_pi, double *OutData_pr, \
			double *OutData_pi, double *out_fading_pr, double *out_fading_pi,\
			int Nr,int Nt,double aglsprR,double *aglR,double aglsprT,double *aglT,\
			double dR,double dT,double *distypeR,double *distypeT,int corrmodel,double alf,\
			double fd,int loson,int uplink)
{ 
  	
    int M_tones = NumOfFreq + 1; 
    long i,tap,Delay_Length;
    double Path_Ampli; 
    double *fading_pr, *fading_pi;
	double *inner_sub_fading_pr,*inner_sub_fading_pi;//
	double *out_sub_fading_pr,*out_sub_fading_pi;//
    long max_delay= (long) Path_Delay[NumOfTaps-1];
    int n;
	int nt,nr;
	double u;
	double pow1,pow2;
	double *hr_r;//
	double *hr_i;
	double *Output_Each_Channel_pr;
	double *Output_Each_Channel_pi;
	double *Output_Sub_path_pr;
	double *Output_Sub_path_pi;
	double *Phase_per_Tap_channel;
	double *matrix_r;
	double *matrix_i;
	double *Hlos_r; //Nr*Nt,Rice H
	double *Hlos_i;
	static bool InData_Complex = 1;//暂时设定为1，待搞清楚

	
    if(corrmodel!=0)
    {
    	fading_pr = (double*)calloc(8,LengthOfBurst);
    	assert(fading_pr);
    	fading_pi = (double*)calloc(8,LengthOfBurst);
    	assert(fading_pi);
		inner_sub_fading_pr = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(inner_sub_fading_pr);
		inner_sub_fading_pi = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(inner_sub_fading_pi);
		out_sub_fading_pr = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(out_sub_fading_pr);
		out_sub_fading_pi = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(out_sub_fading_pi);
		hr_r = (double*)calloc(8,Nr*Nt*Nr*Nt);
		assert(hr_r);
		hr_i = (double*)calloc(8,Nr*Nt*Nr*Nt);
		assert(hr_i);
    }
    else
    {
		fading_pr = (double*)calloc(8,LengthOfBurst);
    	assert(fading_pr);
    	fading_pi = (double*)calloc(8,LengthOfBurst);
    	assert(fading_pi);
		inner_sub_fading_pr = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(inner_sub_fading_pr);
		inner_sub_fading_pi = (double*)calloc(8,Nr*Nt*LengthOfBurst);
		assert(inner_sub_fading_pi);
		out_sub_fading_pr=inner_sub_fading_pr;
		out_sub_fading_pi=inner_sub_fading_pi;
		//hr_r=calloc(8,Nr*Nt*Nr*Nt);
		//assert(hr_r);
		//hr_i=calloc(8,Nr*Nt*Nr*Nt);
		//assert(hr_i);
    }
	if (loson==1)
	{
		matrix_r = (double*)calloc(8,Nr);
		assert(matrix_r);
		matrix_i = (double*)calloc(8,Nr);
		assert(matrix_i);
		Hlos_r = (double*)calloc(8,Nr*Nt);
		assert(Hlos_r);
		Hlos_i = (double*)calloc(8,Nr*Nt);
		assert(Hlos_i);
	}
    
    //
    Output_Each_Channel_pr = (double*)calloc(8,Nr*LengthOfBurst);
    assert(Output_Each_Channel_pr);
    Output_Each_Channel_pi = (double*)calloc(8,Nr*LengthOfBurst);
    assert(Output_Each_Channel_pi);
    Output_Sub_path_pr = (double*)calloc(8,Nr*LengthOfBurst);//
    assert(Output_Sub_path_pr);
    Output_Sub_path_pi = (double*)calloc(8,Nr*LengthOfBurst);
    assert(Output_Sub_path_pi);

    //
    Phase_per_Tap_channel = (double*)calloc(8,M_tones);
    assert(Phase_per_Tap_channel);   
    
	for(i=0;i<max_delay;i++)
	{
		//for(nt=0;nt<Nt;nt++)
	    //{
	    	for(nr=0;nr<Nr;nr++)
	    	{	     		
				Output_Each_Channel_pr[nr+i*Nr] = fore_data_pr[nr+i*Nr];
				Output_Each_Channel_pi[nr+i*Nr] = fore_data_pi[nr+i*Nr];
	        }
	  	//}
   	}
	
	for(i=0;i<Nr*max_delay;i++)
	{
		fore_data_pr[i] = 0.0;
		fore_data_pi[i] = 0.0;
	}
	
	for(nr=0;nr<Nr;nr++)
	{
		for(i=max_delay;i<LengthOfBurst;i++)
		{
			Output_Each_Channel_pr[nr+Nr*i] = 0.0;
			Output_Each_Channel_pi[nr+Nr*i] = 0.0;
		}
	}//clear
	for(i=0;i<Nr*LengthOfBurst;i++)
    {
		Output_Sub_path_pr[i] =0;
		Output_Sub_path_pi[i] =0;
   	}//clear	
    
	for(tap=0;tap<NumOfTaps;tap++)
    {
		for(i=0;i<Nr*LengthOfBurst;i++)
        {
    		Output_Sub_path_pr[i] =0;
    		Output_Sub_path_pi[i] =0;
      	}//clear
		if (loson==1&&tap==0)
		{
			Path_Ampli=Path_Average_Amp[tap];	
			Delay_Length = (long) Path_Delay[tap];//后加，否则没有初始化，导致报错，需要核查是否在los情况下，使用这条语句
			u=2*pi*dT*sin(2*pi*aglT[tap]/360);
     		for(nr=0;nr<Nr;nr++)
    		{
     			matrix_r[nr]=cos(u*nr);
	    		matrix_i[nr]=-sin(u*nr);//
    		}
     		u=2*pi*fd*SampleIndex*dt*cos(alf*2*pi/360);
     		for (nr=0;nr<Nr;nr++)
    		{
	    		for (nt=0;nt<Nt;nt++)
	    		{
	    			if (nt==0)
	    			{
	    				Hlos_r[nr*Nt+nt]=(cos(u)*matrix_r[nr]-sin(u)*matrix_i[nr])*Path_Ampli;
			    		Hlos_i[nr*Nt+nt]=(sin(u)*matrix_r[nr]+cos(u)*matrix_i[nr])*Path_Ampli;
	    			}
		    		else
		    		{
		    			Hlos_r[nr*Nt+nt]=Hlos_r[nr*Nt];
		    			Hlos_i[nr*Nt+nt]=Hlos_i[nr*Nt];
		    		}
  	     		}
	    	}
			for (nr=0;nr<Nr;nr++)
    		{
	    		for (nt=0;nt<Nt;nt++)
	    		{
					for(n=0;n<LengthOfBurst;n++)
	            	{
						out_sub_fading_pr[(nr*Nt+nt)+n*(Nr*Nt)]=Hlos_r[nr*Nt+nt];
						out_sub_fading_pi[(nr*Nt+nt)+n*(Nr*Nt)]=Hlos_i[nr*Nt+nt];
					}
				}
			}
		}
		else
		{
     		for(nr=0;nr<Nr;nr++)
         	{
    	    	for(nt=0;nt<Nt;nt++)
	         	{			    		
	    			Path_Ampli=Path_Average_Amp[tap];	
             	    Delay_Length = (long) Path_Delay[tap];         		        			
             		// get Phase Rotation
        	    	for(i=0;i<M_tones;i++)
	             	{
	            		//Phase_per_Tap_channel[i] = RndPhase[((nr*Nr+nt)*NumOfTaps+tap) +i*NumOfTaps*Nr*Nt];
						Phase_per_Tap_channel[i] = RndPhase[((nr*Nt+nt)*NumOfTaps+tap) +i*NumOfTaps*Nr*Nt];
	             	}                     
    
             		//Generate inner Fading weight
             		JakesIV(Path_Ampli, An, Bn, Wn, Phase_per_Tap_channel, NumOfFreq, \
             			UpdateInterval,UpdatesPerBurst,LengthOfBurst, \
             			(long)SampleIndex, dt, fading_pr, fading_pi);	
					////测试fading第一项是否太大,
					//for (n=0;n<LengthOfBurst;n++)
					//{
					//	if (fading_pi[n] > 1) fading_pi[n] = 0.5;
					//}
    	    	
	             	for(n=0;n<LengthOfBurst;n++)
	            	{
		    			//inner_sub_fading_pr[(nr*Nr+nt)+n*(Nr*Nt)]=fading_pr[n];
		    			//inner_sub_fading_pi[(nr*Nr+nt)+n*(Nr*Nt)]=fading_pi[n];	
						inner_sub_fading_pr[(nr*Nt+nt)+n*(Nr*Nt)]=fading_pr[n];
		    			inner_sub_fading_pi[(nr*Nt+nt)+n*(Nr*Nt)]=fading_pi[n];
		    		}			         												
	    		}
	    	}
		}
		
		if(corrmodel!=0)//
		{
			if (loson==1&&tap==0)
			{
			}
			else
			{
				mimochannel_cx(Nr,Nt,aglsprR,aglR[tap],aglsprT,aglT[tap],dR,dT,(int)distypeR[tap],(int)distypeT[tap],hr_r,hr_i,uplink);
	     		mul_cx(hr_r,hr_i,inner_sub_fading_pr,inner_sub_fading_pi,Nr*Nt,Nr*Nt,LengthOfBurst,out_sub_fading_pr,out_sub_fading_pi);
			}
		}
        
        //generate the output fading
		for(nr=0;nr<Nr;nr++)
     	{
	    	for(nt=0;nt<Nt;nt++)
	     	{
				for(i=0; i<LengthOfBurst; i++)
	            {
					//out_fading_pr[((nr*Nr+nt)*NumOfTaps+tap)+i*NumOfTaps*Nr*Nt]=out_sub_fading_pr[(nr*Nr+nt)+i*Nr*Nt];
					//out_fading_pi[((nr*Nr+nt)*NumOfTaps+tap)+i*NumOfTaps*Nr*Nt]=out_sub_fading_pi[(nr*Nr+nt)+i*Nr*Nt];
					out_fading_pr[((nr*Nt+nt)*NumOfTaps+tap)+i*NumOfTaps*Nr*Nt]=out_sub_fading_pr[(nr*Nt+nt)+i*Nr*Nt];
					out_fading_pi[((nr*Nt+nt)*NumOfTaps+tap)+i*NumOfTaps*Nr*Nt]=out_sub_fading_pi[(nr*Nt+nt)+i*Nr*Nt];
				}
			}
		}
		
		// get the output data of each independent sub-channel
        if(InData_Complex)
    	{
			for(i=0; i<LengthOfBurst; i++)
			{
     			for(nr=0;nr<Nr;nr++)
     			{
     				for(nt=0;nt<Nt;nt++)
     				{
    					//Output_Sub_path_pr[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pr[(nr*Nr+nt)+i*Nr*Nt] 
						//	- InData_pi[nt+i*Nt] * out_sub_fading_pi[(nr*Nr+nt)+i*Nr*Nt]; 
    	            	//Output_Sub_path_pi[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pi[(nr*Nr+nt)+i*Nr*Nt] 
						//	+ InData_pi[nt+i*Nt] * out_sub_fading_pr[(nr*Nr+nt)+i*Nr*Nt]; 
						Output_Sub_path_pr[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pr[(nr*Nt+nt)+i*Nr*Nt] 
							- InData_pi[nt+i*Nt] * out_sub_fading_pi[(nr*Nt+nt)+i*Nr*Nt]; 
    	            	Output_Sub_path_pi[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pi[(nr*Nt+nt)+i*Nr*Nt] 
							+ InData_pi[nt+i*Nt] * out_sub_fading_pr[(nr*Nt+nt)+i*Nr*Nt];
	                }
     			}
    		}
     	}
        else
     	{
    		for(i=0; i<LengthOfBurst; i++)
    		{
    			for(nr=0;nr<Nr;nr++)
    			{
     				for(nt=0;nt<Nt;nt++)
     				{
    					//Output_Sub_path_pr[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pr[(nr*Nr+nt)+i*Nr*Nt]; 
    	            	//Output_Sub_path_pi[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pi[(nr*Nr+nt)+i*Nr*Nt]; 
						Output_Sub_path_pr[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pr[(nr*Nt+nt)+i*Nr*Nt]; 
    	            	Output_Sub_path_pi[nr + i*Nr] += InData_pr[nt+i*Nt] * out_sub_fading_pi[(nr*Nt+nt)+i*Nr*Nt]; 
	                }
		    	}
	    	}
	    }
		
    	for(nr=0;nr<Nr;nr++)
    	{
     		for(i=Delay_Length; i<LengthOfBurst; i++)
            {
     			// adding each sub path to generate outdata
     			Output_Each_Channel_pr[nr + i*Nr] += Output_Sub_path_pr[nr+(i-Delay_Length)*Nr];
    		    Output_Each_Channel_pi[nr + i*Nr] += Output_Sub_path_pi[nr+(i-Delay_Length)*Nr];
    		}
    	}
        
        //conceive the fore frame data
     	//
    	for(nr=0;nr<Nr;nr++)
    	{
    		for(i=0;i<Delay_Length;i++)
    	    {
    			fore_data_pr[nr+i*Nr] += Output_Sub_path_pr[nr+(LengthOfBurst-Delay_Length+i)*Nr];
    		    fore_data_pi[nr+i*Nr] += Output_Sub_path_pi[nr+(LengthOfBurst-Delay_Length+i)*Nr];
    		}
    	}
					
	}//tap
		
	for(nr=0;nr<Nr;nr++)
    {
		for(i=0; i<LengthOfBurst; i++)
    	{
			OutData_pr[nr + i*Nr] = Output_Each_Channel_pr[nr + i*Nr];
    		OutData_pi[nr + i*Nr] = Output_Each_Channel_pi[nr + i*Nr];
	    }
    }
    
	//
	if(Phase_per_Tap_channel) 
		free(Phase_per_Tap_channel);
	if(Output_Sub_path_pr) 
		free(Output_Sub_path_pr);
	if(Output_Sub_path_pi) 
		free(Output_Sub_path_pi);
	if(Output_Each_Channel_pr) 
		free(Output_Each_Channel_pr);
	if(Output_Each_Channel_pi) 
		free(Output_Each_Channel_pi);
	/*if(hr_r) 
		free(hr_r);	
	if(hr_i) 
		free(hr_i);	*/
	if(fading_pr) 
		free(fading_pr);
	if(fading_pi) 
		free(fading_pi);
	if(inner_sub_fading_pr) 
		free(inner_sub_fading_pr);
	if(inner_sub_fading_pi) 
		free(inner_sub_fading_pi);
	if(corrmodel!=0)
	{
		if(out_sub_fading_pr)
			free(out_sub_fading_pr);
		if(out_sub_fading_pi)
			free(out_sub_fading_pi);
		if(hr_r)
			free(hr_r);	
		if(hr_i)
			free(hr_i);
	}
	if (loson==1)
	{
		if (matrix_r)
    	{
    		free(matrix_r);
    	}
    	if (matrix_i)
    	{
    		free(matrix_i);
    	}
     	if(Hlos_r)
     		free(Hlos_r);
    	if(Hlos_i)
    		free(Hlos_i);
	}
	}
	vector<int> ProtocolE::PN9(int length, int seed)
	{
	// 函数名: PN9
	// 函数功能描述:信源数据类型PN9
	// 输入参数: length、seed
	// 输出参数:返回src
	// Called By: GenerateSource
	// 补充说明:
	// 修改日期: 2017/08/08：08:41:41
		// 确保随机种子的范围正确
		if (seed>511)
			seed = 511;
		else if (seed<0)
			seed =0;
		vector<int> src;

		int connection [] ={1,0,0,0,0,1,0,0,0,1};
		int M = 9;// sizeof(connection) / sizeof(connection[0])-1; //寄存器位数
		int con[9]={0};
		for(int i=0;i<M;i++)
			con[i]=connection[i+1]; //用于计算的抽头位数 c1~c9

		int registers[9]={0}; //把10进制种子seed转成二进制比特，填充到寄存器中，低位填充到DM位
		int j=M-1;
		while(seed)
		{
			registers[j]=seed%2;
			seed/=2;
			j--;
		}

		for(int i=0;i<length;i++)
		{
			int tmp=0;
			src.push_back(registers[M-1]); //寄存器M位的值作为PN9序列的一个比特
			for(int m=0;m<M;m++)
				tmp += registers[m] * con[m];//寄存器与抽头系数按位乘

			tmp = tmp % 2;  //模2操作
			for(int n=M-1;n>0;n--)
				registers[n]=registers[n-1]; // 寄存器移位
			registers[0]=tmp; //将tmp2作为反馈的值，输入到寄存器的第一位i
		}

		return src;
	}

	vector<int> ProtocolE:: PN15(int length, int seed)
	{
	// 函数名: PN15
	// 函数功能描述:信源数据类型PN15
	// 输入参数: length、seed
	// 输出参数:返回src
	// Called By: GenerateSource
	// 补充说明:
	// 修改日期: 2017/08/08：08:41:41
	// 输入保护
		if (seed > 32767)
			seed = 32767;
		else if (seed<0)
			seed = 0;

		vector<int> src;

		int connection [] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1};//c14=1，c15=1
		int M = 15;// sizeof(connection) / sizeof(connection[0])-1; //寄存器位数

		int con[15]={0};
		for(int i=0;i<M;i++)
			con[i]=connection[i+1]; //用于计算的抽头位数

		int registers[15]={0};
		int j = M-1;
		while(seed)
		{
			registers[j]=seed%2;
			seed/=2;
			j--;
		}

		for(int i = 0; i < length; i++)
		{
			int tmp=0;
			src.push_back(registers[M-1]); //寄存器M位的值作为PN9序列的一个比特
			for(int m=0;m<M;m++)
				tmp+=registers[m]*con[m];    //寄存器与抽头系数按位乘
			tmp=tmp%2;  //模2操作
			for(int n=M-1;n>0;n--)
				registers[n]=registers[n-1]; // 寄存器移位
			registers[0]=tmp; //将tmp2作为反馈的值，输入到寄存器的第一位i
		}
		return src;
	}

	vector<int> ProtocolE:: All0(int length)
	{
	// 函数名: All0
	// 函数功能描述:信源数据类型PN15
	// 输入参数: length
	// 输出参数:返回src
	// Called By: GenerateSource
	// 补充说明:
	// 修改日期: 2017/08/08：08:43:56
		vector<int> src;
		for (int i = 0; i < length; i++)
			src.push_back(0);
		return src;
	}

	vector<int> ProtocolE:: FromFile(int length, CString filename)
	{
		vector<int> src;
		char* FileName = new char[filename.GetLength()];

		strcpy(FileName,filename.GetBuffer(filename.GetLength()));
		ifstream file;
		file.open(FileName);
		if (!file)
			return src;

		int count = 0;//统计文件位数
		char c = '0';
		while (file.get(c) && count<length) 
		{
			int temp = (c=='0') ? 0 : 1;
			src.push_back(temp);
			count++;
		}
		file.close();

		// 不足位全部补零
		if (count < length)
			for (int i = count; i < length; i++)
				src.push_back(0);

		return src;
	}

	vector<int> ProtocolE:: GenerateSource(int mode,int length, int seed, CString filename)
	{
	// 函数名: GenerateSource
	// 函数功能描述:产生信源数据
	// 输入参数: mode、length、seed、filename
	// 输出参数:返回src
	// Called By: run
	// 补充说明:
	// 修改日期: 2017/08/08：08:47:22
		vector<int> src;
		switch (mode)
		{
		case 0: //PN9
			src = PN9(length,seed);
			break;
		case 1: //PN15
			src = PN15(length,seed);
			break;
		case 2: //ALL0
			src = All0(length);
			break;
		case 3: //read file
			src = FromFile(length,filename);
			break;
		default: //PN9
			src = PN9(length,seed);
			break;
		}
		return src;
	}


	void ProtocolE:: OverSample_Zero(vector<complex<double>>  &Data_In,int Oversampling_Ratio)
	{
	// 函数名: OverSample_Zero
	// 函数功能描述:过采样补零
	// 输入参数: Data_In、Oversampling_Ratio
	// 输出参数:Data_In
	// Called By: run
	// 补充说明:
	// 修改日期: 2017/08/08：08:38:35
		vector<int> zero;
		vector<complex<double>>  Data_temp = Data_In;
		Data_In.clear();
		for (int i=0; i<Oversampling_Ratio-1; i++) zero.push_back(0);
		for (int i=0; i<Data_temp.size(); i++) 
		{
			Data_In.push_back(Data_temp[i]);
			Data_In.insert(Data_In.end(),zero.begin(),zero.end());
		}
		zero.clear();
		Data_temp.clear();
	}

	vector<double> ProtocolE:: Lowpass_filter(int length,double Wn)
	{
	// 函数名: Lowpass_filter
	// 函数功能描述:低通滤波器
	// 输入参数: length、Wn
	// 输出参数:返回filter_coefs
	// Called By: 
	// 补充说明:
	// 修改日期: 2017/08/08：08:46:32
		//win = 0.54-0.46*cos(2*(0:N-1)*pi/(N-1));
		vector<double> filter_coefs;
		vector<double> win;
		vector<double> hd;
		filter_coefs.resize(length);
		hd.resize(length);
		win.resize(length);
		win[0] =0.08;
		for(int i=1;i<length;i++)
			win[i]  =0.54-0.46*cos(2*i*pi/(double)(length-1));
		for(int i=1;i<=length;i++){
			if((i-(length-1)/2.0)!=0)
				hd[i-1]=sin((i-(length-1)/2.0)*Wn*pi)/((i-(length-1)/2.0)*pi);
			else
				hd[i-1] = Wn;
		}
		for(int i=0;i<length;i++){
			filter_coefs[i]=hd[i]*win[i];
		}
		return filter_coefs;
	}
	vector<double> ProtocolE:: Gauss_filter(int length,double BT)
	{
	// 函数名: Gauss_filter
	// 函数功能描述:高斯滤波器
	// 输入参数: length、BT
	// 输出参数:返回filter_coefs
	// Called By: 
	// 补充说明:
	// 修改日期: 2017/08/08：08:45:41
		double Wn = BT;
		vector<double> filter_coefs;
		vector<double> filter_temp;
		vector<double> win;
		vector<double> hd;
		filter_coefs.resize(length);
		hd.resize(length);
		win.resize(length);
		double tao,alpha=3.0,win_max=0.0;
		tao =2.0/Wn;
		for(int i=0;i<length;i++)
			filter_temp.push_back(-2*tao+i*4*tao/(length-1));
		//高斯win=(1/(alpha*sqrt(pi)*tao)).*exp(-(x/(alpha*tao)).^2);
		for(int i=0;i<length;i++)
			win[i]  =(1/(alpha*sqrt(pi)*tao))*exp(-1*pow(filter_temp[i]/(alpha*tao),2.0));
		for(int i=0;i<length;i++)
			if(win[i]>win_max)
				win_max =win[i];
		for(int i=0;i<length;i++)
			win[i]/=win_max;
		//加窗
		for(int i=1;i<=length;i++){
			if((i-(length-1)/2.0)!=0)
				hd[i-1]=sin((i-(length-1)/2.0)*Wn*pi)/((i-(length-1)/2.0)*pi);
			else
				hd[i-1] = Wn;
		}
		for(int i=0;i<length;i++){
			filter_coefs[i]=hd[i]*win[i];
		}
		return filter_coefs;
	}	
	void ProtocolE:: My_conv( vector<complex<double>>  &Data_in,vector<double> &filter_coef,int Data_in_len){
		int len,len_1;
		vector<complex<double>> Data_out;
		vector<complex<double>> Data_out_1;
		complex<double> data_temp;
		double          filter_temp; 
		len =Data_in_len+filter_coef.size()-1;
		len_1 =Data_in_len;
		Data_out.resize(len);
		for (int n=0;n<filter_coef.size();n++)
		{
			for (int i=0;i<=n;i++)
			{
				data_temp =Data_in[n-i];
				filter_temp =filter_coef[i];
				Data_out[n] +=data_temp * filter_temp;
			}
		}

		for (int n=filter_coef.size();n<Data_in_len;n++)
		{
			for (int i=0;i<filter_coef.size();i++)
			{
				data_temp =Data_in[n-i];
				filter_temp =filter_coef[i];
				Data_out[n] +=data_temp * filter_temp;
			}
		}

		for (int n=Data_in_len;n<len+1;n++)
		{
			for (int i=n-Data_in_len+1;i<filter_coef.size();i++)
			{
				data_temp =Data_in[n-i];
				filter_temp =filter_coef[i];
				Data_out[n] +=data_temp * filter_temp;
			}
		}
		for(int i =0 ;i < len_1 ;i++)
		{
			Data_in[i] =Data_out[i];
		}
	}
	vector<double> ProtocolE:: Root_raised_cosine(double alpha,double spb,int ntaps)
	{
	// 函数名: Root_raised_cosine
	// 函数功能描述:根升余弦函数
	// 输入参数: alpha、spb
	// 输出参数:返回taps
	// Called By: run
	// 补充说明:
	// 修改日期: 2017/08/08：08:37:00
		double gain =1.0;
		ntaps |= 1;	// ensure that ntaps is odd
		//double spb = 2; // samples per bit/symbol
		vector<double> taps(ntaps);
		double scale = 0;
		for(int i=0;i<ntaps;i++)
		{
			double x1,x2,x3,num,den;
			double xindx = i - ntaps/2;
			x1 = pi * xindx/spb;
			x2 = 4 * alpha * xindx / spb;
			x3 = x2*x2 - 1;
			if( fabs(x3) >= 0.000001 )  // Avoid Rounding errors...
			{
				if( i != ntaps/2 )
					num = cos((1+alpha)*x1) + sin((1-alpha)*x1)/(4*alpha*xindx/spb);
				else
					num = cos((1+alpha)*x1) + (1-alpha) * pi / (4*alpha);
				den = x3 * pi;
			}
			else
			{
				if(alpha==1)
				{
					taps[i] = -1;
					continue;
				}
				x3 = (1-alpha)*x1;
				x2 = (1+alpha)*x1;
				num = (sin(x2)*(1+alpha)*pi
					- cos(x3)*((1-alpha)*pi*spb)/(4*alpha*xindx)
					+ sin(x3)*spb*spb/(4*alpha*xindx*xindx));
				den = -32 * pi * alpha * alpha * xindx/spb;
			}
			taps[i] = 4 * alpha * num / den;
			scale += taps[i];
		}

		for(int i=0;i<ntaps;i++)
			taps[i] = taps[i] * gain / scale;

		return taps;
	}
	void ProtocolE:: FreOffset_Add( vector<complex<double>>  &Data_in,float FrequencyOffset,float Sampling_Time)
	{
	// 函数名: FreOffset_Add
	// 函数功能描述:频率偏移
	// 输入参数: Data_in、FrequencyOffset、Sampling_Time
	// 输出参数:Data_in
	// Called By: run
	// 补充说明:
	// 修改日期: 2017/08/08：08:39:37
			vector<double> t;
			double Data_temp_real,Data_temp_imag;
			int i;
			for (i=0; i<Data_in.size(); i++ )
			{
				t.push_back(i*Sampling_Time);
			}  

			for (i=0; i<Data_in.size(); i++ )
			{
				Data_temp_real = cos(2*pi*FrequencyOffset*t[i]);
				Data_temp_imag = sin(2*pi*FrequencyOffset*t[i]);
				Data_in[i].real(Data_in[i].real()*Data_temp_real - Data_in[i].imag()*Data_temp_imag);
				Data_in[i].imag(Data_in[i].real()*Data_temp_imag + Data_in[i].imag()*Data_temp_real);
			}  

	}
	//mexErrMsgTxt("right here\n");
