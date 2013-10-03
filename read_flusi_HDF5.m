function [field, box, time] = read_flusi_HDF5( filename )
  %% [field,box] = read_flusi_HDF5( filename )
  % reads in a flusi HDF5 file. the file has to obey standard
  % flusi naming convention. 
  % ----------
  % Input:
  %     filename: the file to be read
  % Output:
  %     field: the actual data
  %     box: the domain size [xl,yl,zl]
  % ----------
  
  
  % open file in read-only mode
  file_id = H5F.open( filename, 'H5F_ACC_RDWR','H5P_DEFAULT');
  
  % open the dataset in the file (in a file mask_00150.h5
  % the dataset is called mask according to our convention)
  dataset_id = H5D.open( file_id,...
               filename(1:strfind(filename,'_')-1)); 
  
  % fetch the resolution from the dataset
  attr_id = H5A.open( dataset_id, 'nxyz');
  nxyz = H5A.read( attr_id );
  H5A.close(attr_id);
  
  % fetch the resolution from the dataset
  attr_id = H5A.open( dataset_id, 'domain_size');
  box = H5A.read( attr_id );
  H5A.close(attr_id);
  
  % fetch time from the dataset
  attr_id = H5A.open( dataset_id, 'time');
  time = H5A.read( attr_id );
  H5A.close(attr_id);
  
  nx = nxyz(1);  
  ny = nxyz(2);
  nz = nxyz(3);  
  fprintf('Resolution is %i %i %i\n',nx,ny,nz)
  
  % read file
  field = H5D.read( dataset_id,'H5ML_DEFAULT', ...  ) % format in memory hier
            'H5S_ALL',... MEMORY    
            'H5S_ALL',... PLATTE
            'H5P_DEFAULT');
  
  
  % close remaining HDF5 objects
  H5D.close(dataset_id)
  H5F.close(file_id)
end