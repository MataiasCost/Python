import os
import pandas as pd
import numpy as np
from os import makedirs
import mplhep as hep
import matplotlib.pyplot as plt
import anatools.data as data
import anatools.analysis as ana
from sklearn.metrics import confusion_matrix
import sys
ana.start()
import gc
import warnings
warnings.filterwarnings('ignore')
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties import umath
from uncertainties.umath import *
import uncertainties
import matplotlib.gridspec as gs
from IPython import get_ipython
import statsmodels.stats.proportion as prop
import math




class HepTrigger:
    def __init__(self,datasets,simulation_samples,data_samples,main_triggers,helper_triggers,channels,variables,label_variaveis,variables_for_SF,regions,data_lumi):
        #Adquirindo as informações dos inputs
        self._datasets = datasets
        self._simulation_samples = simulation_samples
        self._data_samples = data_samples
        self._main_triggers = main_triggers
        self._helper_triggers = helper_triggers
        self._channels = channels
        self._variables = variables
        self._labels = label_variaveis
        self._variables_for_SF = variables_for_SF
        self._dic_info_of_number_of_events_in_each_bin_of_all_triggers = {}
        self._years = self._get_year_samples()
        self._cut_regions = regions
        self._data_lumi = data_lumi
        # Antes de fazer qualquer tipo de estudo, eu vou fazer uma boa pratica quando se esta usando muitos dados
        # Espera-se que todos os nomes dos triggers na lista main_triggers e helper_triggers sejam do tipo bool, ou seja que ocupe apenas 1 byte de memória.
        # Por isso agora eu vou criar uma função voltada para corrigir esse problema de otimização.
        self._change_type_varibles_in_dataframes()


        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue
            #Adquirindo os nomes de todas as samples de dados e simulação
            names_simulation_samples = self._get_datasets_names_of_simulation(year)
            names_data_samples =  self._get_datasets_names_of_data(year)

            #Juntando as amostras e simulações
            datasets[year]['All_data'] = pd.concat( [datasets[year][sample] for sample in names_data_samples]).reset_index(drop=True)
            datasets[year]['All_simulation'] = pd.concat( [datasets[year][sample] for sample in names_simulation_samples]).reset_index(drop=True)
            # Deleting the datasets of MC
            # del datasets['year']['TTTo2L2Nu']
            
            # print(names_data_samples)
            self._list_all_samples = names_simulation_samples + names_data_samples + ['All_data','All_simulation']
            self._list_all_data = names_data_samples + ['All_data']
            self._list_all_simulation = names_simulation_samples + ['All_simulation']
            for sample in self._list_all_samples:
#===========================================================================================================================
# Aplicando correções unicas do minha análise, ou seja, esse bloco de código não é necessário para sua análise! Mas
#caso queira aplicar cortes em variaveis dos seus datasets esse um bom lugar para fazer isso!
# Primeira correção: Aplicar modulo nos valores de duas variaveis: Lep_leading_eta e Lep_subleading_eta
# Segunda correção: Criar as variavéis 'Muon_pt' e 'Electron_pt', que serão utilizadas no canal ElMu
            #Primeira correção
                datasets[year][sample]['Lep_leading_eta'] = datasets[year][sample]['Lep_leading_eta'].abs()
                datasets[year][sample]['Lep_subleading_eta'] = datasets[year][sample]['Lep_subleading_eta'].abs()
            #Segunda correção
                datasets[year][sample]['Muon_pt']=np.where((datasets[year][sample]["RecoLepID"]==1311),datasets[year][sample]["Lep_leading_pt"], (datasets[year][sample]["Lep_subleading_pt"]))
                datasets[year][sample]['Electron_pt']=np.where((datasets[year][sample]["RecoLepID"]==1113),datasets[year][sample]["Lep_leading_pt"], (datasets[year][sample]["Lep_subleading_pt"]))
                datasets[year][sample]['Muon_eta']=np.where((datasets[year][sample]["RecoLepID"]==1311),datasets[year][sample]["Lep_leading_eta"], (datasets[year][sample]["Lep_subleading_eta"]))
                datasets[year][sample]['Electron_eta']=np.where((datasets[year][sample]["RecoLepID"]==1113),datasets[year][sample]["Lep_leading_eta"], (datasets[year][sample]["Lep_subleading_eta"]))

#=============================================================================================================================
            # Vou criar uma nova variável nos datasets reponsável por indicar se o evento ativou algum dos triggers auxiliares ou algum dos triggers principais,
            #essas variáveis vão ser usadas para fazer os graficos e estudos.
            # O ideal é o usuario já implementar isso na lista dos triggers, uma vez que pode ser uma mistura complexa de varias etapas.
            # NO MEU CASO, as variáveis Met_triggers e Lep_triggers já indicam se algum trigger auxiliar ou principal foi ativado, respectivamente. 
                datasets[year][sample]['fired_any_helper_trigger'] = datasets[year][sample]['Met_triggers']
                datasets[year][sample]['fired_any_main_trigger'] = datasets[year][sample]['Lep_triggers']
                # datasets[year][sample]['fired_any_triggers'] = datasets[year][sample]['Met_triggers'] | datasets[year][sample]['Lep_triggers']
                #Criação dos diretorios contendo os graficos (Posso melhorar isso, mas por hora deixo como está):
                # self._create_directorys(year)

                # self.get_all_the_alfas_plots_of_all_triggers(datasets,year)
                                

# Declaração das funções da classe: NÂO PRECISA ALTERAR NADA AQUI que tenh _ antes do nome da função! O resto, se preferir pode alterar (não aconselho, mas eu sou só um aviso também).



    def get_all_SF_systemic_error_of_2D_varibles_for_simulation(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                # gc.collect()
                check_list = []
                for var_x in self._variables_for_SF[channel]:
                    for var_y in self._variables_for_SF[channel]:
                        if var_x == var_y:
                            continue
                        if ( ((str(var_x)+str(var_y)) in check_list) or ((str(var_y)+str(var_x)) in check_list) ):
                            continue
                        check_list.append(str(var_x)+str(var_y))

                        # O erro systematica (até o presente momento) é calculado através de três variaveis:
                        # 1- Diferença entre o Alfa adquirido e o valor teorico dele (valor 1)
                        # Sua ideia veio desse arquiv: https://gitlab.cern.ch/sewuchte/tth-triggerefficiency-dl/-/blob/master/outputs_and_plotting_script/Compute_correlations.C
                        # Essa etapa é facilmente feita usando uma função que já existe:
                        err_Alpha = abs( 1 - self.get_alfa_value('All_simulation',channel,year,'fired_any_main_trigger','fired_any_helper_trigger',show=False))

                        # 2 - Diferença entre o valor do SF e o valor do SF adquirido através da soma das eras:
                        err_eras = self.get_the_diff_between_SF_default_and_SF_of_all_eras(year=year,channel=channel,var_x=var_x,var_y=var_y)

                        # 3 - Maior diferença encontrada em todos os bins para as diferentes regiões:
                        err_regions = self.get_the_diff_between_SF_default_and_SF_of_different_regions(year=year,channel=channel,var_x=var_x,var_y=var_y)


                        print('err_Alpha',err_Alpha)
                        print('err_eras',err_eras)
                        print('err_regions',err_regions)
                        
                        err_systematic = np.zeros(err_eras.shape)
                        for idx_line in range(err_systematic.shape[0]):
                            for idx_col in range(err_systematic.shape[1]):
                                if np.isnan(err_regions[idx_line][idx_col]) or np.isnan(err_eras[idx_line][idx_col]):
                                    continue
                                err_systematic[idx_line][idx_col] = math.sqrt(err_Alpha*err_Alpha +  err_eras[idx_line][idx_col]*err_eras[idx_line][idx_col] + err_regions[idx_line][idx_col]*err_regions[idx_line][idx_col])

                        print(err_systematic)



    def get_the_diff_between_SF_default_and_SF_of_different_regions(self,year,channel,var_x,var_y):
        SF_default = self.get_one_SF_of_2D_varibles(year=year,sample='All_data',channel=channel,var_x=var_x,var_y=var_y,show=False)
        matriz_err = np.zeros(SF_default.shape)
        dict = {}
        for region in self._cut_regions:
            # Adquirindo o valor do SF para cada região
            dict[region] = self.get_one_SF_of_2D_varibles(year=year,sample='All_data',channel=channel,var_x=var_x,var_y=var_y,show=False,region=region)
            # Calculando a diferença de todos os valores, para adquirir a maior diferença
            
            for idx_line in range(dict[region].shape[0]):
                for idx_col in range(dict[region].shape[1]):
                    # print('coordenadas [',idx_line,'][',idx_col,']:',H[idx_line][idx_col])
                    #checa de o valor é igual a NaN (é uma possibilidade)
                    if np.isnan(dict[region][idx_line][idx_col]) or np.isnan(SF_default[idx_line][idx_col]):
                        continue
                    #caso não então a matriz_err pega o maior diferença encontrada entre eles
                    elif abs(SF_default[idx_line][idx_col] - dict[region][idx_line][idx_col]) > matriz_err[idx_line][idx_col]:
                        matriz_err[idx_line][idx_col] = abs(SF_default[idx_line][idx_col] - dict[region][idx_line][idx_col])

        return matriz_err
        

    def get_the_diff_between_SF_default_and_SF_of_all_eras(self,year,channel,var_x,var_y,region='default'):
        SF_default = self.get_one_SF_of_2D_varibles(year=year,sample='All_data',channel=channel,var_x=var_x,var_y=var_y,show=False,region=region)
        dict = {}
        for sample in self._data_lumi[year]:
            dict[sample] = self.get_one_SF_of_2D_varibles(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,show=False,region=region)
        
        #Agora é aplicado um fator em cada era
        for key,value in dict.items():
            # print('value:',value)
            # print('É aplicado:',self._data_lumi[year][key])
            dict[key] = value*self._data_lumi[year][key]
            # print('ficou:',dict[key])

        #Somando as eras em apenas uma variavel
        SF_eras = np.sum( SF_era for SF_era in dict.values())
        return abs(SF_eras-SF_default)

    def get_one_SF_estatistic_error_of_2D_varibles_for_simulation(self,year,channel,var_x,var_y,type,show=False):


        data_var_x_bin_changed = self._variables[channel][var_x][:]
        data_var_x_bin_changed[-1] = self._datasets[year]['All_data'][var_x].max()
        data_var_y_bin_changed = self._variables[channel][var_y][:]
        data_var_y_bin_changed[-1] = self._datasets[year]['All_data'][var_y].max()

        simu_var_x_bin_changed = self._variables[channel][var_x][:]
        simu_var_x_bin_changed[-1] = self._datasets[year]['All_data'][var_x].max()
        simu_var_y_bin_changed = self._variables[channel][var_y][:]
        simu_var_y_bin_changed[-1] = self._datasets[year]['All_data'][var_y].max()

        # self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True')[var_x].reset_index(drop=True)
        Data_numerator, _, _,_ = plt.hist2d(
                                            x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[data_var_x_bin_changed,data_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )
        # Data_selected_2, _, _,_ = plt.hist2d(
        #                                     x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
        #                                     ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
        #                                     ,bins=[self._variables[channel][var_x],self._variables[channel][var_y]]
        #                                     ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
        #                                     )
        Data_denominator, _, _,_ = plt.hist2d(
                                        x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[data_var_x_bin_changed,data_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )


        # Data_not_selected_2, _, _,_ = plt.hist2d(
        #                                     x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
        #                                     ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
        #                                     ,bins=[self._variables[channel][var_x],self._variables[channel][var_y]]
        #                                     ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
        #                                     )




    #samples[1] = MC_el_met
        MC_selected, _, _,_ = plt.hist2d(
                                            x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )


        MC_selected_2, _, _,_ = plt.hist2d(
                                            x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )




        MC_not_selected, _, _,_ = plt.hist2d(
                                            x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )



        MC_not_selected_2, _, _,_ = plt.hist2d(
                                            x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
                                            ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
                                            ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                            ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                            )



        #Adquirimos o valor do r e calculamos o seu erro
        # r_data = Data_not_selected/Data_selected
        # r_err_data = r_data*np.sqrt( (Data_not_selected_2/Data_not_selected)**2  + (Data_selected_2/Data_selected)**2 ) 
        # Calculando o erro da SIMULAÇÂO
        r_mc = MC_not_selected/MC_selected
        r_err_mc = r_mc*np.sqrt( (MC_not_selected_2/MC_not_selected)**2  + (MC_selected_2/MC_selected)**2 ) 


        simulation_eff_ratio = 1/(1+r_mc)
        #propagação de erro
        simulation_err = (r_err_mc)/np.power(1+r_mc,2)

        # print('simulation_err',simulation_err)
        # print('r_err_mc',r_err_mc)

        data_eff_err = []
        for x in range(r_err_mc.shape[1]):
            for y in range(r_err_mc.shape[0]):
                data_eff_err.append(r_err_mc[y,x]/( 1+r_mc[y,x] )**2 ) 

        # Calculando o erro dos DADOS
        # https://www.statsmodels.org/stable/generated/statsmodels.stats.proportion.proportion_confint.html
        y_below, y_above = prop.proportion_confint(Data_numerator, Data_denominator, alpha=0.32, method='beta')
        data_eff_ratio = np.divide(Data_numerator,Data_denominator)
        data_err_below = data_eff_ratio - y_below
        data_err_above = y_above - data_eff_ratio
        # print('ratio',ratio)
        # print('Err_static_eff_MC_'+channel+'_'+var_x+'_'+var_y+':',ye_below ,' e ', ye_above)

        # print('data_eff_ratio',data_eff_ratio)
        # print('data_err_above',data_err_above)

        if type == 'upper':
            error = np.sqrt( np.power(simulation_err/simulation_eff_ratio,2) + np.power(data_err_above/data_eff_ratio,2)  )
        if type == 'bottom':
            error = np.sqrt( np.power(simulation_err/simulation_eff_ratio,2) + np.power(data_err_below/data_eff_ratio,2)  )
        
        # print('error',error)

        hist=np.transpose(error)
        plt.close()

        fig, ax = plt.subplots()
        ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
        p = ax.pcolormesh(hist)
            
        for i in range(len(hist)):
            for j in range(len(hist[i])):
                text = ax.text(j+0.5, i+0.5, round(hist[i, j], 3), color="w", ha="center", va="center", fontweight="bold")
            
            
        ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
        ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])

        ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
        ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
        
        ax.set_xlabel(self._labels[var_x])
        ax.set_ylabel(self._labels[var_y])

        # plt.close()
        
        if show == False:
            plt.close()
        
        return error












    def do_all_SF_estatistic_error_of_2D_varibles_for_simulation(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                # gc.collect()
                types = ['upper','bottom']
                for type in types:
                    check_list = []
                    for var_x in self._variables_for_SF[channel]:
                        for var_y in self._variables_for_SF[channel]:
                            if var_x == var_y:
                                continue
                            if ( ((str(var_x)+str(var_y)) in check_list) or ((str(var_y)+str(var_x)) in check_list) ):
                                continue
                            check_list.append(str(var_x)+str(var_y))


                            #Ideia burra que eu aprendi a fazer, deve ter uma forma melhor de fazer isso

                            data_var_x_bin_changed = self._variables[channel][var_x][:]
                            data_var_x_bin_changed[-1] = self._datasets[year]['All_data'][var_x].max()
                            data_var_y_bin_changed = self._variables[channel][var_y][:]
                            data_var_y_bin_changed[-1] = self._datasets[year]['All_data'][var_y].max()

                            simu_var_x_bin_changed = self._variables[channel][var_x][:]
                            simu_var_x_bin_changed[-1] = self._datasets[year]['All_data'][var_x].max()
                            simu_var_y_bin_changed = self._variables[channel][var_y][:]
                            simu_var_y_bin_changed[-1] = self._datasets[year]['All_data'][var_y].max()
    ###################################
                            # self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True')[var_x].reset_index(drop=True)
                            Data_numerator, _, _,_ = plt.hist2d(
                                                                x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[data_var_x_bin_changed,data_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )


                            # Data_selected_2, _, _,_ = plt.hist2d(
                            #                                     x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                            #                                     ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                            #                                     ,bins=[self._variables[channel][var_x],self._variables[channel][var_y]]
                            #                                     ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                            #                                     )




                            Data_denominator, _, _,_ = plt.hist2d(
                                                            x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[data_var_x_bin_changed,data_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )


                            # Data_not_selected_2, _, _,_ = plt.hist2d(
                            #                                     x = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
                            #                                     ,y = self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
                            #                                     ,bins=[self._variables[channel][var_x],self._variables[channel][var_y]]
                            #                                     ,weights=self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_data'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                            #                                     )




                        #samples[1] = MC_el_met
                            MC_selected, _, _,_ = plt.hist2d(
                                                                x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )


                            MC_selected_2, _, _,_ = plt.hist2d(
                                                                x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == True').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )




                            MC_not_selected, _, _,_ = plt.hist2d(
                                                                x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )



                            MC_not_selected_2, _, _,_ = plt.hist2d(
                                                                x = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_x].reset_index(drop=True)
                                                                ,y = self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])[var_y].reset_index(drop=True)
                                                                ,bins=[simu_var_x_bin_changed,simu_var_y_bin_changed]
                                                                ,weights=self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)*self._datasets[year]['All_simulation'].query('fired_any_helper_trigger == True and fired_any_main_trigger == False').query(self._channels[channel])['evtWeight'].reset_index(drop=True)
                                                                )



                            #Adquirimos o valor do r e calculamos o seu erro
                            # r_data = Data_not_selected/Data_selected
                            # r_err_data = r_data*np.sqrt( (Data_not_selected_2/Data_not_selected)**2  + (Data_selected_2/Data_selected)**2 ) 
                            # Calculando o erro da SIMULAÇÂO
                            r_mc = MC_not_selected/MC_selected
                            r_err_mc = r_mc*np.sqrt( (MC_not_selected_2/MC_not_selected)**2  + (MC_selected_2/MC_selected)**2 ) 


                            simulation_eff_ratio = 1/(1+r_mc)
                            #propagação de erro
                            simulation_err = (r_err_mc)/np.power(1+r_mc,2)

                            # print('simulation_err',simulation_err)
                            # print('r_err_mc',r_err_mc)

                            data_eff_err = []
                            for x in range(r_err_mc.shape[1]):
                                for y in range(r_err_mc.shape[0]):
                                    data_eff_err.append(r_err_mc[y,x]/( 1+r_mc[y,x] )**2 ) 

                            # Calculando o erro dos DADOS
                            # https://www.statsmodels.org/stable/generated/statsmodels.stats.proportion.proportion_confint.html
                            y_below, y_above = prop.proportion_confint(Data_numerator, Data_denominator, alpha=0.32, method='beta')
                            data_eff_ratio = np.divide(Data_numerator,Data_denominator)
                            data_err_below = data_eff_ratio - y_below
                            data_err_above = y_above - data_eff_ratio
                            # print('ratio',ratio)
                            # print('Err_static_eff_MC_'+channel+'_'+var_x+'_'+var_y+':',ye_below ,' e ', ye_above)

                            # print('data_eff_ratio',data_eff_ratio)
                            # print('data_err_above',data_err_above)

                            if type == 'upper':
                                error = np.sqrt( np.power(simulation_err/simulation_eff_ratio,2) + np.power(data_err_above/data_eff_ratio,2)  )
                            if type == 'bottom':
                                error = np.sqrt( np.power(simulation_err/simulation_eff_ratio,2) + np.power(data_err_below/data_eff_ratio,2)  )
                            
                            # print('error',error)

                            hist=np.transpose(error)
                            plt.close()

                            fig, ax = plt.subplots()
                            ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                            p = ax.pcolormesh(hist)
                                
                            for i in range(len(hist)):
                                for j in range(len(hist[i])):
                                    text = ax.text(j+0.5, i+0.5, round(hist[i, j], 3), color="w", ha="center", va="center", fontweight="bold")
                                
                                
                            ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
                            ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
        
                            ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
                            ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                            
                            ax.set_xlabel(self._labels[var_x])
                            ax.set_ylabel(self._labels[var_y])

                            # plt.close()
                            
                            
                            os.makedirs('./'+year+'/2D/Err_SF_statistic',exist_ok=True)
                            plt.savefig('./'+year+'/2D/Err_SF_statistic/Err_static_'+type+'_of_eff_'+channel+'_'+var_x+'_'+var_y+'.png') 
                            plt.close()
                            





    def get_one_SF_of_2D_varibles(self,year,sample,channel,var_x,var_y,show=False,region='default'):

        Eff_data = self.get_one_eff_of_2D_varibles(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,show=False,region=region)
        Eff_simulation = self.get_one_eff_of_2D_varibles(year=year,sample='All_simulation',channel=channel,var_x=var_x,var_y=var_y,show=False,region=region)
        
        
        
        ratio = np.divide(Eff_data,Eff_simulation)
        fig, ax = plt.subplots()
        ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
        hist=np.transpose(ratio)
        p = ax.pcolormesh(hist)
        
        for i in range(len(hist)):
            for j in range(len(hist[i])):
                text = ax.text(j+0.5, i+0.5, round(hist[i, j], 4), color="w", ha="center", va="center", fontweight="bold")
                
        ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
        ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
        
        ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
        ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                            
        ax.set_xlabel(self._labels[var_x])
        ax.set_ylabel(self._labels[var_y])

        if show == False:
            plt.close()

        del Eff_data,Eff_simulation,p,text,ax,fig
        gc.collect()
        return ratio



    def do_all_SF_of_2D_varibles(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_data:
                # gc.collect()
                    check_list = []
                    for var_x in self._variables_for_SF[channel]:
                        for var_y in self._variables_for_SF[channel]:
                            if var_x == var_y:
                                continue
                            if ( ((str(var_x)+str(var_y)) in check_list) or ((str(var_y)+str(var_x)) in check_list) ):
                                continue
                            check_list.append(str(var_x)+str(var_y))

                            # print('uepaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')

                            # print(bins)
                            # print('Estamos indo no year:',year,'sample:',sample,'channel',channel)
                            Eff_data = self.get_one_eff_of_2D_varibles(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,show=False)
                            Eff_simulation = self.get_one_eff_of_2D_varibles(year=year,sample='All_simulation',channel=channel,var_x=var_x,var_y=var_y,show=False)
                            
                            
                            
                            # print(eff)

                            # print('Os bins utilizados foram: ',bins_x_metXlep, ' e ',bins_y_metXlep)
                            
                            ratio = np.divide(Eff_data,Eff_simulation)

                            fig, ax = plt.subplots()
                            ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                            hist=np.transpose(ratio)
                            p = ax.pcolormesh(hist)
                            
                            for i in range(len(hist)):
                                for j in range(len(hist[i])):
                                    text = ax.text(j+0.5, i+0.5, round(hist[i, j], 4), color="w", ha="center", va="center", fontweight="bold")
                            
                            # ax_x = 
                            
                            ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
                            ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
                            
                            ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
                            ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                                                
                            ax.set_xlabel(self._labels[var_x])
                            ax.set_ylabel(self._labels[var_y])

                            
                            # print(ratio)
                            # plt.savefig('ueapa.png')
                            os.makedirs('./'+year+'/2D/SF_of_bins',exist_ok=True)
                            plt.savefig('./'+year+'/2D/SF_of_bins/SF_'+sample+'_'+channel+'_'+var_x+'_'+var_y+'.png') 
                            plt.close()

                            del ratio,Eff_data,Eff_simulation,p,text,ax,fig

                    gc.collect()
                            # return ratio,bins_x_met,bins_y_met




    def get_one_eff_of_2D_varibles(self,year,sample,channel,var_x,var_y,show=False,region='default'):
        N_metXlep = self.get_2D_number_of_events_per_bin(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,events_that_already_fired_what_trigg='helperAndmain',region=region)
        N_met = self.get_2D_number_of_events_per_bin(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,events_that_already_fired_what_trigg='helper',region=region)
        plt.close()
        # eff = []

        # for idx in range(len(N_metXlep)):
        #     eff.append(N_metXlep[idx]/N_met[idx])

        # print(eff)

        # print('Os bins utilizados foram: ',bins_x_metXlep, ' e ',bins_y_metXlep)

        ratio = np.divide(N_metXlep,N_met)

        fig, ax = plt.subplots()
        ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
        hist=np.transpose(ratio)
        p = ax.pcolormesh(hist)

        for i in range(len(hist)):
            for j in range(len(hist[i])):
                text = ax.text(j+0.5, i+0.5, round(hist[i, j], 4), color="w", ha="center", va="center", fontweight="bold")

        # ax_x = 

        ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
        ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])

        ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
        ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                    
        ax.set_xlabel(self._labels[var_x])
        ax.set_ylabel(self._labels[var_y])


        # print(ratio)
        # plt.savefig('ueapa.png')
        # os.makedirs('./'+year+'/2D/Eff_of_bins',exist_ok=True)
        # plt.savefig('./'+year+'/2D/Eff_of_bins/Eff_'+sample+'_'+channel+'_'+var_x+'_'+var_y+'.png') 
        if show == False:
            plt.close()

        return ratio
    
    
    #ESSA FUNÇÃO NÃO FUNCIONA AINDA, ESTOU TENTANDO FORMULAR ELA MELHOR
    def my_reset(self, *varnames):
        """
        varnames são as variáveis que você deseja manter, por default utilizo ele nos atributos e funções das classes
        """
        globals_ = globals()
        to_save = {v: globals_[v] for v in varnames}
        functions = [f for f in dir(self) if callable(getattr(self, f))]
        for f in functions:
            to_save[f] = getattr(self, f)
        attributes = [a for a in dir(self) if not callable(getattr(self, a))]
        for a in attributes:
            to_save[a] = getattr(self, a)
        to_save['my_reset'] = self.my_reset
        to_save['objeto'] = self
        del globals_
        get_ipython().magic('reset -sf') 
        
        globals().update(to_save)



    def _change_type_varibles_in_dataframes(self):
        #Primeiro vou criar uma lista de todos os termos que devem ser salvos, os demais vão ser excluidos
        for year in self._years:
            if year != '2018':
                continue
            for sample in (self._data_samples[year] + self._simulation_samples[year]):
                # print(sample,year)
            # mudando o tipo de variavel dos main_triggers:
                for helper_trigger in self._helper_triggers[year]:
                   self._datasets[year][sample][helper_trigger] = self._datasets[year][sample][helper_trigger].astype(bool)
                for main_trigger in self._main_triggers[year]:
                   self._datasets[year][sample][main_trigger] = self._datasets[year][sample][main_trigger].astype(bool)



    def do_all_eff_of_2D_varibles(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_samples:
                    # gc.collect()
                    check_list = []
                    for var_x in self._variables_for_SF[channel]:
                        for var_y in self._variables_for_SF[channel]:
                            if var_x == var_y:
                                continue
                            if ( ((str(var_x)+str(var_y)) in check_list) or ((str(var_y)+str(var_x)) in check_list) ):
                                continue
                            check_list.append(str(var_x)+str(var_y))

                            # print('uepaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')

                            # print(bins)
                            # print('Estamos indo no year:',year,'sample:',sample,'channel',channel)
                            N_metXlep = self.get_2D_number_of_events_per_bin(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,events_that_already_fired_what_trigg='helperAndmain')
                            N_met = self.get_2D_number_of_events_per_bin(year=year,sample=sample,channel=channel,var_x=var_x,var_y=var_y,events_that_already_fired_what_trigg='helper')
                            
                            # eff = []
                            
                            # for idx in range(len(N_metXlep)):
                            #     eff.append(N_metXlep[idx]/N_met[idx])
                            
                            # print(eff)

                            # print('Os bins utilizados foram: ',bins_x_metXlep, ' e ',bins_y_metXlep)
                            
                            ratio = np.divide(N_metXlep,N_met)

                            fig, ax = plt.subplots()
                            ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                            hist=np.transpose(ratio)
                            p = ax.pcolormesh(hist)
                            
                            for i in range(len(hist)):
                                for j in range(len(hist[i])):
                                    text = ax.text(j+0.5, i+0.5, round(hist[i, j], 4), color="w", ha="center", va="center", fontweight="bold")
                            
                            # ax_x = 
                            
                            ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
                            ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
                            
                            ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
                            ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                                                
                            ax.set_xlabel(self._labels[var_x])
                            ax.set_ylabel(self._labels[var_y])

                            
                            # print(ratio)
                            # plt.savefig('ueapa.png')
                            os.makedirs('./'+year+'/2D/Eff_of_bins',exist_ok=True)
                            plt.savefig('./'+year+'/2D/Eff_of_bins/Eff_'+sample+'_'+channel+'_'+var_x+'_'+var_y+'.png') 
                            plt.close()
                            # return ratio,bins_x_met,bins_y_met
    
    

    #Teste Unitario, será implementado de forma correta mais pra frente do framework
    # Resultado: VALIDADO
    def _TEST_get_2D_number_of_events_per_bin(self):
        values_trigg = [0,0,1,0,1,0,0,1]
        values_lead = [20,40,45,70,90,120,150,30]
        values_sublead = [10,30,40,69,40,110,120,15]
        x_bin = [0,30,60,90,120,150]
        y_bin = [0,60,120,150]

        x = values_trigg[:]
        x_var = values_lead[:]
        y = values_trigg[:]
        y_var = values_sublead[:]

        for indice in range(len(x)):
            if x[indice] == True:
                x[indice] = x_var[indice]
            elif x[indice] == False:
                x[indice] = -999
        for indice in range(len(y)):
            if y[indice] == True:
                y[indice] = y_var[indice]
            elif y[indice] == False:
                y[indice] = -999
        
        hist, _, _, _ = plt.hist2d(x,y,bins=[x_bin,y_bin])
        print(hist)
        plt.close()
        list = []            
        for i in range(len(y_bin)-1):
            for j in range(len(x_bin)-1):
                list.append(int(hist.T[i,j]))
        
        fig, ax = plt.subplots()
        ratio=np.transpose(hist)
                            # print(ratio)
        p = ax.pcolormesh(ratio)
                            
        for i in range(len(ratio)):
            for j in range(len(ratio[i])):
                text = ax.text(j+0.5, i+0.5, round(ratio[i, j], 2), color="w", ha="center", va="center", fontweight="bold")

        ax.set_yticks(np.arange(0,len(y_bin)))
        ax.set_yticklabels([str(bin) for bin in y_bin])
                            
        ax.set_xticks(np.arange(0,len(x_bin)))
        ax.set_xticklabels([str(bin) for bin in x_bin])

        # return hist


    def get_2D_number_of_events_per_bin(self,year,sample,channel,var_x,var_y,show=False,events_that_already_fired_what_trigg='helperAndmain',region='default'):
        plt.figure(figsize=(16,12))
        cut_region = self._cut_regions[region]
        #########################################################################################################
        #Devido ao fato de algumas regiões estarem atreladas ao canal que elas se encontram, é necessário implementar ela aqui:
        z=0
        if region == 'nMuon0' and channel == 'MuMu':
            z=2
        else:
            z=0
        if region == 'nMuon1' and channel == 'MuMu':
            z=3
        else:
            z=1
        if region == 'nMuon2' and channel == 'MuMu':
            z=4
        else:
            z=2
        if region == 'nElectron0' and channel == 'ElEl':
            z=2
        else:
            z=0
        if region == 'nElectron1' and channel == 'ElEl':
            z=3
        else:
            z=1
        if region == 'nElectron2' and channel == 'ElEl':
            z=4
        else:
            z=2    
        if channel== 'ElMu':
            z=0

        #########################################################################################################
        # print('CHEGOUUUUUUUUUUUUUU')
        ax2 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
        ana.style(ax2, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
        if events_that_already_fired_what_trigg == 'helperAndmain':
            # print('CHEGOUUUUUUUUUUUUUU21')
            x = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True').query(cut_region.format(z))['fired_any_main_trigger'].tolist()
            # print('CHEGOUUUUUUUUUUUUUU22')
            x_var = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True').query(cut_region.format(z))[var_x].tolist()
            # print('CHEGOUUUUUUUUUUUUUU23')
            y = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True').query(cut_region.format(z))['fired_any_main_trigger'].tolist()
            # print('CHEGOUUUUUUUUUUUUUU24')
            y_var = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True').query(cut_region.format(z))[var_y].tolist()
        elif events_that_already_fired_what_trigg == 'helper':
            x = self._datasets[year][sample].query(self._channels[channel]).query(cut_region.format(z))['fired_any_helper_trigger'].tolist()
            x_var = self._datasets[year][sample].query(self._channels[channel]).query(cut_region.format(z))[var_x].tolist()
            y = self._datasets[year][sample].query(self._channels[channel]).query(cut_region.format(z))['fired_any_helper_trigger'].tolist()
            y_var = self._datasets[year][sample].query(self._channels[channel]).query(cut_region.format(z))[var_y].tolist()
        else:
            print('Essa função é feita para calcular a eficiência nos graficos 2D, por isso só é possível ter events_that_already_fired_what_trigg \n '+ \
                  +'= helper, o que vai calcular o número de eventos que passaram pelo helepr e pelo main ou = \'no one\' que representa o número de eventos que passaram só pelo helper trigg')
        # print('CHEGOUUUUUUUUUUUUUU2')

        for indice in range(len(x)):
            if x[indice] == True:
                x[indice] = x_var[indice]
            elif x[indice] == False:
                x[indice] = -999
        for indice in range(len(y)):
            if y[indice] == True:
                y[indice] = y_var[indice]
            elif y[indice] == False:
                y[indice] = -999        

        var_x_bins_changed = self._variables[channel][var_x][:]
        var_y_bins_changed = self._variables[channel][var_y][:]
        var_x_bins_changed[-1] = self._datasets[year][sample][var_x].max()
        var_y_bins_changed[-1] = self._datasets[year][sample][var_y].max()


        if events_that_already_fired_what_trigg == 'helperAndmain':
            hist, _, _, _ = plt.hist2d(x,y,
                                       bins=[var_x_bins_changed,var_y_bins_changed],
                                       weights=(self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True').query(cut_region.format(z))['evtWeight'].tolist()))
        if events_that_already_fired_what_trigg == 'helper':
            hist, _, _, _ = plt.hist2d(x,y,
                                       bins=[var_x_bins_changed,var_y_bins_changed],
                                       weights=(self._datasets[year][sample].query(self._channels[channel]).query(cut_region.format(z))['evtWeight'].tolist()))
        
        # print('CHEGOUUUUUUUUUUUUUU3')
        plt.xticks(self._variables[channel][var_x], [str(bin) for bin in self._variables[channel][var_x]])
        ax2.set_xlabel(self._labels[var_x])
        plt.yticks(self._variables[channel][var_y], [str(bin) for bin in self._variables[channel][var_y]],rotation=90)
        ax2.set_ylabel(self._labels[var_y])
        plt.colorbar()
        plt.close()

        # print('CHEGOUUUUUUUUUUUUUU4')
        list = []            
        for i in range(len(self._variables[channel][var_y])-1):
            for j in range(len(self._variables[channel][var_x])-1):
                list.append(int(hist.T[i,j]))
                                
        # plt.close()
                            # # print(int(hist.T[0,0]))
                            
        _, ax = plt.subplots()
        ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                            
                            # print(hist)
        ratio=np.transpose(hist)
                            # print(ratio)
        p = ax.pcolormesh(ratio)
                            
        for i in range(len(ratio)):
            for j in range(len(ratio[i])):
                ax.text(j+0.5, i+0.5, round(ratio[i, j], 2), color="w", ha="center", va="center", fontweight="bold")
                            
        ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
        ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
                            
        ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
        ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                            
        ax.set_xlabel(self._labels[var_x])
        ax.set_ylabel(self._labels[var_y])

        if show == False:                
            plt.close() 


        #Vou fazer um teste, e deletar todas as variáveis que usei nessa função:
        del ratio,p,x_var,x,y,y_var
        # print('Limpando...')
        gc.collect()




        return hist

####################################################################################################################################
    #FALTA ACRESCENTAR O ULTIMO BIN = np.inf!
    def do_2D_of_all_number_of_events_per_bin(self):
        #Para impedir que haja duplicatas de algum histograma, mas com os eixos trocados, eu vou implementar uma lista que vai decidir se as duas variaveis ja foram usadas
        # Ela vai se chamar check_list
        list_what_trigg_the_data_already_fired=['helperAndmain','helper']
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_samples:
                    for context in list_what_trigg_the_data_already_fired:
                        check_list = []
                        # print(context)
                        for var_x in self._variables_for_SF[channel]:
                            for var_y in self._variables_for_SF[channel]:
                                if var_x == var_y:
                                    continue
                                if ( ((str(var_x)+str(var_y)) in check_list) or ((str(var_y)+str(var_x)) in check_list) ):
                                    continue
                                check_list.append(str(var_x)+str(var_y))
                                fig1 = plt.figure(figsize=(16,12))

                                ax2 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
                                ana.style(ax2, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                                
                                
                                if context == 'helperAndmain':
                                    x = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True')['fired_any_main_trigger'].tolist()
                                    x_var = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True')[var_x].tolist()
                                    y = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True')['fired_any_main_trigger'].tolist()
                                    y_var = self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True')[var_y].tolist()
                                elif context == 'helper':
                                    x = self._datasets[year][sample].query(self._channels[channel])['fired_any_helper_trigger'].tolist()
                                    x_var = self._datasets[year][sample].query(self._channels[channel])[var_x].tolist()
                                    y = self._datasets[year][sample].query(self._channels[channel])['fired_any_helper_trigger'].tolist()
                                    y_var = self._datasets[year][sample].query(self._channels[channel])[var_y].tolist()
                                else:
                                    print('''Essa função é feita para calcular a eficiência nos graficos 2D, por isso só é possível ter events_that_already_fired_what_trigg \n 
                                        = helper, o que vai calcular o número de eventos que passaram pelo helepr e pelo main ou = \'no one\' que representa o número de eventos que passaram só pelo helper trigg''')
                                    

                                for indice in range(len(x)):
                                    if x[indice] == True:
                                        x[indice] = x_var[indice]
                                    elif x[indice] == False:
                                        x[indice] = -999
                                for indice in range(len(y)):
                                    if y[indice] == True:
                                        y[indice] = y_var[indice]
                                    elif y[indice] == False:
                                        y[indice] = -999        

                                var_x_bins_changed = self._variables[channel][var_x][:]
                                var_y_bins_changed = self._variables[channel][var_y][:]
                                var_x_bins_changed[-1] = self._datasets[year][sample][var_x].max()
                                var_y_bins_changed[-1] = self._datasets[year][sample][var_y].max()
                                
                                if context == 'helperAndmain':
                                    hist, _, _, _ = plt.hist2d(x,y,bins=[var_x_bins_changed,var_y_bins_changed],weights=(self._datasets[year][sample].query(self._channels[channel]).query('fired_any_helper_trigger == True')['evtWeight'].tolist()))
                                if context == 'helper':
                                    hist, _, _, _ = plt.hist2d(x,y,bins=[var_x_bins_changed,var_y_bins_changed],weights=(self._datasets[year][sample].query(self._channels[channel])['evtWeight'].tolist()))
                                                                

                                plt.xticks(self._variables[channel][var_x], [str(bin) for bin in self._variables[channel][var_x]])
                                ax2.set_xlabel(self._labels[var_x])
                                plt.yticks(self._variables[channel][var_y], [str(bin) for bin in self._variables[channel][var_y]],rotation=90)
                                ax2.set_ylabel(self._labels[var_y])
                                plt.colorbar()
                                plt.close()


                                list = []            
                                #================================================================================================================================
                                for i in range(len(self._variables[channel][var_y])-1):
                                    for j in range(len(self._variables[channel][var_x])-1):
                                        # plt.text((bins[var[0]][j+1]+bins[var[0]][j])/2,(bins[var[1]][i+1]+bins[var[1]][i])/2, int(hist.T[i,j]), color="w", ha="center", va="center", fontweight="bold")
                                        list.append(int(hist.T[i,j]))
                                    
                                plt.close()
                                # # print(int(hist.T[0,0]))
                                
                                fig, ax = plt.subplots()
                                ana.style(ax, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                                
                                # print(hist)
                                ratio=np.transpose(hist)
                                # print(ratio)
                                p = ax.pcolormesh(ratio)
                                
                                for i in range(len(ratio)):
                                    for j in range(len(ratio[i])):
                                        text = ax.text(j+0.5, i+0.5, round(ratio[i, j], 2), color="w", ha="center", va="center", fontweight="bold")
                                
                                ax.set_yticks(np.arange(0,len(self._variables[channel][var_y])))
                                ax.set_yticklabels([str(bin) for bin in self._variables[channel][var_y]])
                                
                                ax.set_xticks(np.arange(0,len(self._variables[channel][var_x])))
                                ax.set_xticklabels([str(bin) for bin in self._variables[channel][var_x]])
                                
                                ax.set_xlabel(self._labels[var_x])
                                ax.set_ylabel(self._labels[var_y])

                                os.makedirs('./'+year+'/2D/Number_of_events',exist_ok=True)

                                plt.savefig('./'+year+'/2D/Number_of_events/N_evt_'+sample+'_'+channel+'_'+var_x+'_'+var_y+'_that_pass_'+context+'.png')

                            
                                plt.close() 



    def do_all_eff_trigger_of_simulation_in_comparation_with_data(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_simulation:
                    for var in self._variables[channel]:
                        fig1 = plt.figure(figsize=(8,6))
                        grid = [2, 1]
                        gs1 = gs.GridSpec(grid[0], grid[1], height_ratios=[4, 1])
                        N = 1
                        ax1 = plt.subplot() 
                        ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1

                        # print(channel)
                        # print(var)
                        # print(self._variables[channel][var])
                        #Testando se tem mais de uma sample:
                        yratio_data, _, _ = ana.efficiency_plot( ax1, var, self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel]), "fired_any_main_trigger", label='Data', color='b',histograms=False, bins=self._variables[channel][var],overflow = True)
                        yratio_simulation, _, _ = ana.efficiency_plot( ax1, var, self._datasets[year][sample].query(self._channels[channel]), "fired_any_main_trigger", label='Simulation', color='r',histograms=False,weight="evtWeight", bins=self._variables[channel][var],overflow = True)
                        ana.style(ax1, lumi=59.8, year=year, legend_ncol=1, legend_loc='center right',ylim=[0.6,1.1],xticks=self._variables[channel][var])
                        ana.labels(ax1, ylabel="Events")  # Set up the label names
                        # ana.labels(ax1, xlabel=self._labels[var], ylabel=r'Efficiency')  # Set up the label names
                        # ax1.text(0.35, 0.01, 'Channel '+ channel,verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='green', fontsize=15)
                        

                        # ybkg, errbkg = ana.stacked_plot( ax1, var, dataframes, labels, colors, weight="evtWeight", bins=bins )  # Produce the stacked plot
                        # ydata, errdata = ana.data_plot( ax1, var, df_DATA, bins=bins )
                        
                        #==================================================
                        ax2 = plt.subplot(ana.position(gs1,grid,N,2), sharex=ax1)  # Positioning at subplot 2 of the plot number 2
                        #==================================================
                        ana.ratio_plot( ax2, yratio_data, np.zeros(len(yratio_data)), yratio_simulation, np.zeros(len(yratio_simulation)), bins=self._variables[channel][var])
                        ana.labels(ax2, xlabel=self._labels[var], ylabel="Data / Bkg.")  # Set up the label names
                        ana.style(ax2, ylim=[0.95, 1.05], yticks=[0.95, 1, 1.05], xgrid=True, ygrid=True,xticks=self._variables[channel][var])
                        
                        
                        os.makedirs('./'+year+'/Efficiency',exist_ok=True)
                        plt.savefig('./'+year+'/Efficiency/Eff_'+sample+'_'+channel+'_'+var+'.png')     
                        
                        
                        
                        # plt.savefig('uepa.png')
                        plt.close()
                        # print('uepa')

    def get_eff_trigger_of_simulation_and_data(self,year,sample,channel,var,show=False):
        fig1 = plt.figure(figsize=(8,6))
        grid = [2, 4]
        gs1 = gs.GridSpec(grid[0], grid[1], height_ratios=[1, 1])
        N = 1
        ax1 = plt.subplot()
        yratio_data, _, _ = ana.efficiency_plot( ax1, var, self._datasets[year]['All_data'].query('fired_any_helper_trigger == True').query(self._channels[channel]), "fired_any_main_trigger", label='Data', color='r',histograms=False, bins=self._variables[channel][var],overflow = True)
        yratio_simulation, _, _ = ana.efficiency_plot( ax1, var, self._datasets[year][sample].query(self._channels[channel]), "fired_any_main_trigger", label='Simulation', color='b',histograms=False,weight="evtWeight", bins=self._variables[channel][var],overflow = True)
        ana.style(ax1, lumi=59.8, year=year, legend_ncol=1, legend_loc='center right',ylim=[0.6,1.1])
        ana.labels(ax1, xlabel=self._labels[var], ylabel=r'Efficiency')  # Set up the label names
        ax1.text(0.35, 0.01, 'Channel '+ channel,verticalalignment='bottom', horizontalalignment='right',transform=ax1.transAxes,color='green', fontsize=15)
        if show==False:
            plt.close()
        return yratio_data,yratio_simulation
        
                        
    def do_all_the_plots_of_N_eff(self):
        types = ['pass','not pass']
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_simulation:
                    for var in self._variables[channel]:
                        for type in types:
                            N_eff_after,N_eff_not_after = self._get_N_eff_of_the_bins(sample=sample,channel=channel,year=year,var=var)
                            plt.close()    
                            bins_2 = list()
                            valores = []
                            #Corrigindo os valores NaN
                            # print('Valor do N_eff after: ',N_eff_after)
                            N_eff_after = [x if x==x else 0 for x in N_eff_after]
                            # print('Valor do N_eff after: ',N_eff_after)
                            if (type == "pass"):
                                for indice, elemento in enumerate(N_eff_after):
                                    bins_2.append(indice)
                                    # print(elemento)
                                    valores.append(int(elemento))
                                    tag = "win"
                            N_eff_not_after = [x if x==x else 0 for x in N_eff_not_after]
                            if (type == "not pass"):
                                for indice, elemento in enumerate(N_eff_not_after):
                                    bins_2.append(indice)
                                    # print(elemento)
                                    valores.append(int(elemento))
                                    tag = "lose"
                                
                            fig,ax = plt.subplots()#plt.figure(figsize=(8,6))
                            ax.xaxis.tick_bottom()
                            plt.tick_params(
                                axis='x',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom=False,      # ticks along the bottom edge are off
                                top=False,         # ticks along the top edge are off
                                labelbottom=True) # labels along the bottom edge are off
                            plt.tick_params(
                                axis='y',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right=False,         # ticks along the top edge are off
                                labelbottom=True    # labels along the bottom edge are off
                            ) 
                            plt.axis('on')
                            ax.bar(bins_2, valores, width=1, color='red',edgecolor='blue',align= 'edge')
                            ax.yaxis.set_ticklabels([])
                            plt.ylabel(r"$N_{eff}$ "+type)
                            plt.xlabel(self._labels[var])
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            ax.spines['left'].set_visible(False)
                            bins_2.append(bins_2[-1]+1)
                            # print(bins_2)
                            # print(bins)
                            plt.xticks(bins_2, self._variables[channel][var], rotation='horizontal',size='small')
                            ax.tick_params(labeltop = False) 
                            for i, v in enumerate(valores):
                                if (v==0):
                                    plt.text(i,v + 1, 'NaN', color="blue", fontweight='bold')
                                else:
                                    plt.text(i,v + 1, str(v), color="blue", fontweight='bold')
                            # plt.show()
                            #==============================================================================================================================================
                            os.makedirs('./'+year+'/N_eff',exist_ok=True)
                            plt.savefig('./'+year+'/N_eff/N_eff_'+sample+'_'+channel+'_'+var+'_'+tag+'.png')
                            # plt.close()
                    
    def _get_N_eff_of_the_bins(self,year,sample,channel,var):
                        
        datasets_selected = self._datasets[year][sample][self._datasets[year][sample]['fired_any_main_trigger'] == 1].query(self._channels[channel]).query('fired_any_helper_trigger == True')
        datasets_not_selected = self._datasets[year][sample][self._datasets[year][sample]['fired_any_main_trigger'] == 0].query(self._channels[channel]).query('fired_any_helper_trigger == True')

        # overflow and underflow
        eff_bins = self._variables[channel][var][:]
        eff_bins[-1] = np.inf
        # eff_bins[0] = -np.inf
                        
                        
        # y_before, _ = np.histogram( datasets[year][sample][var], eff_bins, weights=datasets[year][sample]['evtWeight'] )
        y_after, _ = np.histogram( datasets_selected[var], eff_bins, weights=datasets_selected['evtWeight'] )    
                        
                        
        y_not_after, _ = np.histogram( datasets_not_selected[var], eff_bins, weights=datasets_not_selected['evtWeight'] ) 
        y2_not_after, _ = np.histogram( datasets_not_selected[var], eff_bins, weights=datasets_not_selected['evtWeight']*datasets_not_selected['evtWeight'] )
        # y2_before, _ = np.histogram( datasets[year][sample][var], eff_bins, weights=datasets[year][sample]['evtWeight']*datasets[year][sample]['evtWeight'] )
        y2_after, _ = np.histogram( datasets_selected[var], eff_bins, weights=datasets_selected['evtWeight']*datasets_selected['evtWeight'] )


        # err_after = np.sqrt(y2_after)
        # err_not_after = np.sqrt(y2_not_after)
                            
        # r = y_not_after/y_after
        # err_r = r*np.sqrt( (err_after/y_after)**2 + (err_not_after/y_not_after)**2 )
        # err_ratio = err_r/(1+r)**2
                            
        N_eff_after = y_after**2/y2_after
        N_eff_not_after = y_not_after**2/y2_not_after

        return N_eff_after,N_eff_not_after 


    def get_alfa_value(self,sample,channel,year,trigg_1,trigg_2,show=False):
        x = self._datasets[year][sample].query(self._channels[channel])[trigg_1].tolist()
        for indice in range(len(x)):
            if x[indice] == True:
                x[indice] = 1
            elif x[indice] == False:
                x[indice] = 0
        y = self._datasets[year][sample].query(self._channels[channel])[trigg_2].tolist()
        for indice in range(len(y)):
            if y[indice] == True:
                y[indice] = 1
            elif y[indice] == False:
                y[indice] = 0
        label_bins_x_and_y = ['Not Fired', 'Fired']
        hist, xbins, ybins, _ = plt.hist2d(x,y,bins=[0,1,2],weights=(self._datasets[year][sample].query(self._channels[channel])['evtWeight'].tolist()))
        # print(hist)

        plt.xticks([0.5,1.5], label_bins_x_and_y)
        plt.text(1, -0.20, trigg_1, ha='center')
        plt.text(0.05, 1.95,sample+', Channel '+channel, va='center', rotation=0,color='r',size=20)
        plt.yticks([0.5,1.5], label_bins_x_and_y,rotation=90)
        
        plt.text(-0.3, 1, trigg_2, va='center', rotation='vertical')
        plt.colorbar()
        # ybins = [0,1,2]
        # xbins = [0,1,2]
    #############################################################################################################
        # plt.text(0.40, 1, r"$\alpha$: "+str(round(alfa4.n,4))+r"$\pm$"+str(round(err_alfa,4)),size=25,color= 'w', va='center', rotation=0)
    #############################################################################################################
        for i in range(len(ybins)-1):
            for j in range(len(xbins)-1):
                plt.text(xbins[j]+0.5,ybins[i]+0.5, round(hist.T[i,j],1), color="w", ha="center", va="center", fontweight="bold")


        # print(hist)
        eff_main = (hist[1][0]+hist[1][1])/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
        eff_helper = (hist[0][1]+hist[1][1])/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
        eff_of_both_main_and_helper = hist[1][1]/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
        value_alfa = (eff_main*eff_helper)/eff_of_both_main_and_helper
        # print(value_alfa)
        plt.text(0.7, 1, r"$\alpha$: "+str(round(value_alfa,4)),size=25,color= 'w', va='center', rotation=0) #)+r"$\pm$"+str(round(err_alfa,4)

        if show == False:
            plt.close()

        
        
        return value_alfa



    def do_all_the_alfas_plots_of_all_triggers(self):
        for year in self._years:
            #Segurança enquanto só tem implementado o ano de 2018
            if year != '2018':
                continue    
            for channel in self._channels:
                for sample in self._list_all_samples:
                    # print('SAMPLE: ',sample)
                    for helper_trigg in (self._helper_triggers[year] + ['fired_any_helper_trigger']):
                        # print(helper_trigg)
                        # print('Sample:',sample,' Channel:', channel)
                        # ax2 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
                        # ana.style(ax2, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                        label_bins_x_and_y = ['Not Fired', 'Fired']
                        x = self._datasets[year][sample].query(self._channels[channel])['fired_any_main_trigger'].tolist()
                        for indice in range(len(x)):
                            if x[indice] == True:
                                x[indice] = 1
                            elif x[indice] == False:
                                x[indice] = 0
                        y = self._datasets[year][sample].query(self._channels[channel])[helper_trigg].tolist()
                        for indice in range(len(y)):
                            if y[indice] == True:
                                y[indice] = 1
                            elif y[indice] == False:
                                y[indice] = 0
                        
                        ax2 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
                        ana.style(ax2, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                        hist, xbins, ybins, im = plt.hist2d(x,y,bins=[0,1,2],weights=(self._datasets[year][sample].query(self._channels[channel])['evtWeight'].tolist()))
                        # print(hist)
                        # print(len(x))
                        # print(len(datasets[data].loc[datasets[data]['Lep_triggers']==True]['evtWeight'].tolist()))
                        # [datasets[data].loc[['Lep_triggers']==True]['evtWeight'],datasets[data].loc[[MET_trigger]==True]['evtWeight']]
                        plt.xticks([0.5,1.5], label_bins_x_and_y)
                        plt.text(1, -0.20, 'Any main trigger', ha='center')
                        plt.text(0.05, 1.95,sample+', Channel '+channel, va='center', rotation=0,color='r',size=20)
                        plt.yticks([0.5,1.5], label_bins_x_and_y,rotation=90)
                        if (helper_trigg == 'fired_any_helper_trigger'):
                            plt.text(-0.3, 1, 'Any helper trigger', va='center', rotation='vertical')
                        else:
                            plt.text(-0.3, 1, helper_trigg, va='center', rotation='vertical')
                        plt.colorbar()
                        # ybins = [0,1,2]
                        # xbins = [0,1,2]
                    #############################################################################################################
                        # plt.text(0.40, 1, r"$\alpha$: "+str(round(alfa4.n,4))+r"$\pm$"+str(round(err_alfa,4)),size=25,color= 'w', va='center', rotation=0)
                    #############################################################################################################
                        for i in range(len(ybins)-1):
                            for j in range(len(xbins)-1):
                                plt.text(xbins[j]+0.5,ybins[i]+0.5, round(hist.T[i,j],1), color="w", ha="center", va="center", fontweight="bold")


                        # print(hist)
                        eff_main = (hist[1][0]+hist[1][1])/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
                        eff_helper = (hist[0][1]+hist[1][1])/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
                        eff_of_both_main_and_helper = hist[1][1]/(hist[0][0] + hist[1][0]+ hist[0][1]+ hist[1][1])
                        value_alfa = (eff_main*eff_helper)/eff_of_both_main_and_helper
                        # print(value_alfa)
                        plt.text(0.7, 1, r"$\alpha$: "+str(round(value_alfa,4)),size=25,color= 'w', va='center', rotation=0) #)+r"$\pm$"+str(round(err_alfa,4)

                        os.makedirs('./'+year+'/Plots_Alfas',exist_ok=True)
                        plt.savefig('./'+year+'/Plots_Alfas/Alfa_'+sample+'_'+channel+'_'+helper_trigg+'.png')
                        plt.close()
                        # return alfa4.n





    def do_hist_of_number_of_events_pass_of_each_helper_trigger(self):
        types = ['values','percent']
        for year in self._years:
            if year != '2018':
                continue    
            for sample in self._list_all_samples:    
                for channel in self._channels.keys():    
                    for type in types: 
                        # hep.cms.label('Work in progress', data=True, lumi=59.8, year=year)
                        ax1 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
                        ana.style(ax1, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                        for idx in range(len(self._helper_triggers[year])):
                            xbins = 0+idx*0.25
                            if type == 'values':
                                plt.bar(xbins,self._datasets[year][sample].query(self._channels[channel]).loc[self._datasets[year][sample][self._helper_triggers[year][idx]]==True]['evtWeight'].sum(),width=0.25,label=self._helper_triggers[year][idx])
                            if type == 'percent':
                                plt.bar(xbins,self._datasets[year][sample].query(self._channels[channel]).loc[self._datasets[year][sample][self._helper_triggers[year][idx]]==True]['evtWeight'].sum()/self._datasets[year][sample].query(self._channels[channel])['evtWeight'].sum(),width=0.25,label=self._helper_triggers[year][idx])
                        plt.legend(loc='upper left',fontsize=15)
                        if type == 'values':
                            plt.ylabel("Firing Rate")
                        if type == 'percent':
                            plt.ylim(0,1)
                            plt.ylabel("Percent")
                        plt.xticks([])
                        plt.xlabel(sample+', channel '+channel)  
                        os.makedirs('./'+year+'/trigger_info/',exist_ok=True)
                        if type == 'values':
                            plt.savefig('./'+year+'/trigger_info/'+'Values_'+sample+'_'+channel+'_helper_trigg.png')
                        if type == 'percent':
                            plt.savefig('./'+year+'/trigger_info/'+'Percent_'+sample+'_'+channel+'_helper_trigg.png')
                        plt.close()


    def do_hist_of_number_of_events_pass_of_each_main_trigger(self):
        types = ['values','percent']
        for year in self._years:
            if year != '2018':
                continue
            for sample in self._list_all_samples:
                for channel in self._channels.keys():
                    for type in types:           
                        # hep.cms.label('Work in progress', data=True, lumi=59.8, year=year)
                        ax1 = plt.subplot()              # Positioning at subplot 1 of the plot number 1
                        ana.style(ax1, lumi=59.8, year=year, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top
                        for idx in range(len(self._main_triggers[year])):
                            xbins = 0+idx*0.25
                            if type == 'values':
                                plt.bar(xbins,self._datasets[year][sample].query(self._channels[channel]).loc[self._datasets[year][sample][self._main_triggers[year][idx]]==True]['evtWeight'].sum(),width=0.25,label=self._main_triggers[year][idx])
                            if type == 'percent':
                                plt.bar(xbins,self._datasets[year][sample].query(self._channels[channel]).loc[self._datasets[year][sample][self._main_triggers[year][idx]]==True]['evtWeight'].sum()/self._datasets[year][sample].query(self._channels[channel])['evtWeight'].sum(),width=0.25,label=self._main_triggers[year][idx])
                        plt.legend(loc='upper left',fontsize=15)
                        if type == 'values':
                            plt.ylabel("Firing Rate")
                        if type == 'percent':
                            plt.ylim(0,1)
                            plt.ylabel("Percent")
                        plt.xticks([])
                        plt.xlabel(sample+', channel '+channel)  
                        os.makedirs('./'+year+'/trigger_info/',exist_ok=True)
                        if type == 'values':
                            plt.savefig('./'+year+'/trigger_info/'+'Values_'+sample+'_'+channel+'_main_trigg.png')
                        if type == 'percent':
                            plt.savefig('./'+year+'/trigger_info/'+'Percent_'+sample+'_'+channel+'_main_trigg.png')
                        plt.close()




    def _create_directorys(self,path):
        os.makedirs('./'+path,exist_ok=True)

    def _get_year_samples(self):
        return self._simulation_samples.keys()
        
    def _get_datasets_names_of_simulation(self,year):
        return self._simulation_samples[year]
    
    def _get_datasets_names_of_data(self,year):
        return self._data_samples[year]


        
        
