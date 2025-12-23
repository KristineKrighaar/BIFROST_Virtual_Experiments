import numpy as np
import h5py
import pandas as pd
import warnings
import os
import pickle

back_info = pd.read_csv('BIFROST_McStas_backend_information.csv')
#back_info = back_info.loc[back_info['wedge_number'] == 1]

back_info = np.asarray(back_info)[:,0:5]

for j in range(len(back_info)):
    back_info[j,3] = back_info[j,3].replace(',','.')
    back_info[j,4] = back_info[j,3].replace(',','.')

back_info = np.asarray(back_info, dtype=float)    
    
back_info = back_info[(back_info[:,0]== 2.7)  & (back_info[:,2]== 4) | (back_info[:,0]== 3.2) & (back_info[:,2]== 4)| (back_info[:,0]== 3.8) & (back_info[:,2]== 5) | (back_info[:,0]== 4.4) & (back_info[:,2]== 5) | (back_info[:,0]== 5.0) & (back_info[:,2]== 5)]


class BIFROST_measurement:
    """
    This Object is designed to contain 1 tube mearsurment, in a given (A3, A4) setting,  including data and metadata needed for further reduction. To this object the following information can be assosiated:

    **Should be defined from:** wedge, arc, tube

    Data:
    - I (y, t)
    - I _err (y,t)

    Axis: t_s, y_m, A4

    Metadata:
    - A3 (sample rotation)
    - L_sd (length from detector to sample)
    - Ef (The Ef for the given analyser)
    """
    def __init__(self,wedge=0, arc=0, tube=0, I=0, I_err=0, t_s = 0, y_m = 0, A4=0, A3=0, L_sd=0, Ef=0, Delta_E=0, Ei=0, dA4=0, Q=0, L_sdi=0):
        self.wedge = wedge
        self.arc = arc
        self.tube = tube
        self.I = I
        self.I_err = I_err
        self.t_s = t_s
        self.y_m = y_m
        self.A3 = A3
        self.A4 = A4
        self.Delta_E = Delta_E
        self.Ei = Ei
        self.dA4 = dA4
        self.Q = Q
        
        self.L_sd = L_sd
        self.Ef = Ef
    
        """ Get L_sd """
        arc_values = {0: 2.7,1: 3.2,2: 3.8,3: 4.4,4: 5.0}

        # Loop through each arc value
        for arc, value in arc_values.items():
            # Check if the current arc value matches self.arc
            if self.arc == arc:
                # Filter the corresponding row from back_info
                ar = back_info[back_info[:, 0] == value]
                
                # Determine wedge conditions and calculate self.L_sd accordingly
                if self.wedge in (0, 3, 6):
                    self.L_sd = ar[0, 3] + ar[0, 4]
                elif self.wedge in (1, 4, 7):
                    self.L_sd = ar[1, 3] + ar[1, 4]
                elif self.wedge in (2, 5, 8):
                    self.L_sd = ar[2, 3] + ar[2, 4]
                break  # Exit the loop after finding the matching arc value

        """ Get Ef (Think about making the Ef correction based on the dA4 takeoff angle) """
        if arc == 0:
            self.Ef = 2.7
        if arc == 1:
            self.Ef = 3.2
        if arc == 2:
            self.Ef = 3.8
        if arc == 3:
            self.Ef = 4.4
        if arc == 4:
            self.Ef = 5.0

            
    def getMcStasData(self, filename):
        """
        Loading all relevant data from McStas simulation mccode file.

        Takes: string with file path to relevant folder

        **returns**
        Attributes for: 
        - I: Intensities in shape (# y bins, # t bins) 
        - I_err: Intensities errors, same shape as I
        - t_s: ToF bins on detector.
        - y_m: Spacial pixels on detector. Given in distance from center. shape (100)
        - A3: integer for sample rotation
        - dA4: Angular offset from center of detectortube to each pixel. Same shpe as y_m
        - A4: Angular offset in scattering angle from direct beam to individual pixel.
        """
    
        I = np.loadtxt(str(filename)+'/signal_1Dspace_'+str(self.wedge)+'_'+str(self.arc)+'_'+str(self.tube)+'.dat').reshape(3,1000,100)[0]
        self.I = I.T
        I_err = np.loadtxt(str(filename)+'/signal_1Dspace_'+str(self.wedge)+'_'+str(self.arc)+'_'+str(self.tube)+'.dat').reshape(3,1000,100)[1]
        self.I_err = I_err.T

        file = str(filename)+'/signal_1Dspace_'+str(self.wedge)+'_'+str(self.arc)+'_'+str(self.tube)+'.dat'

        with open(file, 'r') as the_file:
            all_data = [line.strip() for line in the_file.readlines()]
            limits = all_data[28]
            limits = limits.replace('# xylimits: ', '')
            limits = limits.split(' ')
            limits = np.asarray(limits, dtype=float)
            self.t_s = np.linspace(limits[2], limits[3], 1000)
            self.y_m = np.linspace(limits[0], limits[1], 100)
    
            A4 = all_data[9]
            A4 = A4.replace('# Param: A4=', '')
            A4_offset = np.asarray(A4).astype('float')

            A3 = all_data[10]
            A3 = A3.replace('# Param: A3=', '')
            self.A3 = np.asarray(A3).astype('float') 


        ################## Modify A4 for each pixel in tube #####################
        # Get the angular offset from wedge number
        wedge_offsets = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80])
        wedge_offset = np.array(wedge_offsets[self.wedge])
            
        # Calculate the A4 degrees for each pixel on tube
        self.dA4 = np.degrees(np.arctan((self.y_m/self.L_sd)))

        # Add the offset for tank and wedge number
        A4 = (A4_offset + wedge_offset + self.dA4)*-1
        self.A4 = np.flip(A4)

        return self.I, self.I_err, self.t_s, self.y_m, self.A3, self.A4, self.dA4

    def calcDE(self):
        """
        Calculates the energy transfer for each measured data-point I.
        1) Determines the Ef_offset at for each detectorpixel along the tube.
        2) Calculates Delta_E for each Ef (dA4), T. 
        3) Returns Ef and Ei to have the same shape (len(t), len(y)) as Delta_E since it is needed for calculating Q.

        **return** Delta_E, Ef, Ei in same shape as I.

        """
        
        # Make correction to each physicals bins Ef (At this step self.Ef has shape (100,))
        ideal_lam = 1/(0.11056*np.sqrt(self.Ef))
        self.Ef = (9.045/(ideal_lam*np.cos(np.radians(self.dA4))))**2

        # Define constants
        t_offset_optimal = 6940.56*1e-6
        t = self.t_s - t_offset_optimal # time measured in [s] has shape (N time bins)
        Ef_J = self.Ef*1.602176634e-22 # [J]
        Li = (162-6.35)# [m]
        m_n = 1.67492749804e-27 # [kg]

        # Make the correction to the travel path in the secondary spectrometer
        """
        1) Lave et array som regner hvad længden er til hver pixel på en detektor. L_sd = sin(dA4)*y_m
        2) Make sure that the array vertical and then stack in the horizontal direction 1000 times so we have a 100x1000 matrix.
        """
        self.L_sd = abs(self.y_m)/np.sin(np.radians(abs(self.dA4)))
        L_sd = np.tile(self.L_sd[:, np.newaxis], (1, len(t)))

        # calculate Delta E for all points I
        T, D = np.meshgrid(t, Ef_J)
        #print('T_shape =',np.shape(T))
        #print(D)
        self.Ei = (Li**2*m_n)/(2*(T-np.sqrt((m_n*L_sd**2)/(2*D)))**2)
        self.Delta_E = (self.Ei-D)*6.24150907e21 # convert from J to meV
        
        self.Ei = self.Ei*6.24150907e21 # convert from J to meV

        # Make Ef the right shape to have Ef for each I datapoint
        self.Ef = self.Ef.reshape(-1,1) 
        self.Ef = np.tile(self.Ef, (1,len(t)))

        return self.Delta_E, self.Ef, self.Ei # has shape (y_m=100, N time bins)
    

    def calcQ(self):
        """
        
        """
        # Calc ki
        lam_i = 1/(0.11056*np.sqrt(self.Ei))
        ki = (2*np.pi)/lam_i
        
        # Calc kf
        lam_f = 1/(0.11056*np.sqrt(self.Ef))
        kf =  ((2*np.pi)/lam_f)
        #print('Shape of kf =',np.shape(kf))
        
        # Reshape A4 to have A4 for each matrix
        A4 = self.A4.reshape(-1,1)
        A4 = np.tile(A4, (1,len(ki[0,:])))
    
        # Calc qx and qy       
        qx = ki-kf*np.cos(np.radians(A4))
        qy = -kf*np.sin(np.radians(A4))

        Q = np.stack((qx,qy)).reshape(2,-1)

        # Convert A3 to radians
        A3 = np.radians(self.A3)

        # Define the Sample rotation matrix
        Rot_matrix = np.array([[np.cos(A3), -np.sin(A3)],[np.sin(A3), np.cos(A3)]])

        self.Q = np.dot(Rot_matrix, Q)
        return self.Q

    def flatten(self):
        """
        
        """
        self.I = self.I.flatten()

        self.I_err = self.I_err.flatten()

        self.Delta_E = self.Delta_E.flatten()

        qx = self.Q[0].flatten()
        qy = self.Q[1].flatten()
        self.Q = np.stack((qx,qy))
        return self.I, self.I_err, self.Delta_E, self.Q

    def CorrectI(self):
        """
        
        """
        m_n = 1.67492749804e-27 # [kg]
    
        de_dt = (m_n*(162-6.35+self.L_sd)**2)/(self.t_s**3)*6.24150907e21 # convert from J to meV
        
        # Reshape A4 to have A4 for each matrix
        de_dt = np.tile(de_dt, (len(self.I[:,0]),1))
        self.I = self.I/de_dt
        
        return self.I


    def Normalize(self):
        """
        Should take the tube data and normalize to the integrated intensity in each pixel.
        """
        with open("vana_norm.pkl", "rb") as file:
            loaded_obj = pickle.load(file)
        
        return 1


    def load_Backend_measurment(filepath, Correct=False, Norm_data=False):
        """
        
        """
        setting = {}
        
        I = []
        I_err = []
        DeltaE = []
        qx = []
        qy = []

        
        for w in range(9):
            #for t in range(3):
                obj = BIFROST_measurement(wedge=w, arc=4, tube=1)

                obj.getMcStasData(filepath)
                obj.calcDE()
                obj.calcQ()
                if Correct==True:
                    obj.CorrectI()
                
                if Norm_data==False:
                    obj.flatten()
                    I.append(obj.I)
                    I_err.append(obj.I_err)
                    DeltaE.append(obj.Delta_E)
                    qx.append(obj.Q[0])
                    qy.append(obj.Q[1])

                setting['TM_'+str(w)+'_4_'+str(1)] = obj
        
        if Norm_data==True:
            with open("vana_norm.pkl", "wb") as file:
                pickle.dump(setting, file)

        I = np.asarray(I).flatten()
        I_err = np.asarray(I_err).flatten()
        DeltaE = np.asarray(DeltaE).flatten()
        qx = np.asarray(qx).flatten()
        qy = np.asarray(qy).flatten()

        return I, I_err, DeltaE, qx, qy
    
    def load_scan(keyword, output_filename, folder_path='.', reduce_zero=False, set_Elim=False, mcgui=False, Correct=False, Norm_data=False):
        """
        Load scan data by searching for files containing a keyword in a specified folder.

        Parameters:
            keyword (str): Keyword to look for in filenames.
            output_filename (str): Name of the output file to create or use.
            folder_path (str): Path to the folder to search in. Default is current directory.
            reduce_zero (bool): Optional flag.
            set_Elim (bool): Optional flag.
            mcgui (bool): Optional flag.
            Correct (bool): Optional flag.
            Norm_data (bool): Optional flag.
        """
         # Ensure the folder exists
        if not os.path.isdir(folder_path):
            print(f"Folder '{folder_path}' does not exist.")
            return
        
        # Search for files containing the keyword in the specified folder
        files_with_keyword = [file for file in os.listdir(folder_path) if keyword in file]

        if not files_with_keyword:
            print(f"No files found with the specified keyword '{keyword}' in '{folder_path}'.")
            return
        
        # If simulation is from mcgui, then the files_with_keyword should be found in folders
        if mcgui == True:
            files_with_keyword = [os.path.join(keyword, name) for name in os.listdir(keyword) if os.path.isdir(os.path.join(keyword, name))]

        # 
        with open(output_filename, 'a') as output_file:
            for file_with_keyword in files_with_keyword:
                
                file_with_keyword = folder_path+file_with_keyword

                # Load information from the file
                processed_data = BIFROST_measurement.load_Backend_measurment(file_with_keyword, Correct=Correct, Norm_data=Norm_data)
                
                #print('shape original =',np.shape(processed_data[0]))

                # Convert tuple of results to nd.array
                matrix = np.column_stack(processed_data)


                # Reduce all datapoints where I = 0
                if reduce_zero == True:
                    matrix = matrix[matrix[:,0]>0]
                
                print('Energier = ', matrix[:,2])

                # Choose the upper and lower energy bound saved to txt file.
                if set_Elim != False:
                    if set_Elim[0]>set_Elim[1]:
                        warnings.warn("Lower bound is greater than upper bound. Change order of numbers.", category=Warning)
                    matrix = matrix[(matrix[:,2]>=set_Elim[0]) & (matrix[:,2]<=set_Elim[1])]

                print('shape matrix =',np.shape(matrix))
                # Append processed data to output file
                np.savetxt(output_file, matrix, delimiter='\t', fmt='%f')
                output_file.write('\n')

    