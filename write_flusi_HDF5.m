function write_flusi_HDF5( filename, field, box, time )
  %% write_flusi_HDF5( filename, field, domain )
  % dumps a field in a flusi filename according to the internal
  % convention. The resolution+domain size are stored in attributes. 
  % You may ignore the domain size and pass [0 0 0] instead, depending
  % on what the field is used for.
  % You must store one field per file only.
  % The filename convention is NAME_TIME.h5, e.g. mask_00100.h5
  % ----------
  % Input: 
  %     filename: the file to be written to
  %     field: the actual data
  %     box: the box size of the domain [xl,yl,zl] (you may set [0 0 0])
  %     time: timestamp of the snapshot (you may set 0.0)
  % Output:
  %     none
  % ----------
  
  % dataset name
  dsetname = filename( 1:strfind(filename,'_')-1 );
  
  % create file
  file_id=H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

  % create the dataspace (note C vs FORTRAN array indexing style)
  d = 3; % 3D data
  dataspace_id = H5S.create_simple( d, fliplr(size(field)),...
                                    fliplr(size(field))  );

  % create dataset
  dataset_id = H5D.create( file_id, dsetname, 'H5T_IEEE_F32LE',...
               dataspace_id,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');

  %-------------------------          
  % attributes 
  %-------------------------
  % resolution
  acpl_id  = H5P.create('H5P_ATTRIBUTE_CREATE');
  type_id  = H5T.copy('H5T_STD_I32LE');
  space_id = H5S.create_simple(1,3,3);  
  attr_id  = H5A.create(dataset_id,'nxyz',type_id,space_id,acpl_id);  
  H5A.write ( attr_id, type_id, int32(size(field)) )
  
  % domain size 
  acpl_id  = H5P.create('H5P_ATTRIBUTE_CREATE');
  type_id  = H5T.copy('H5T_IEEE_F32LE');
  space_id = H5S.create_simple(1,3,3);  
  attr_id  = H5A.create(dataset_id,'domain_size',type_id,space_id,acpl_id);  
  H5A.write ( attr_id, type_id, single(box) )  
  
  % timestamp 
  acpl_id  = H5P.create('H5P_ATTRIBUTE_CREATE');
  type_id  = H5T.copy('H5T_IEEE_F32LE');
  space_id = H5S.create_simple(1,1,1);  
  attr_id  = H5A.create(dataset_id,'time',type_id,space_id,acpl_id);  
  H5A.write ( attr_id, type_id, single(time) ) 
  
  %-------------------------
  % dump dataset 
  %-------------------------
  H5D.write ( dataset_id,'H5ML_DEFAULT', ...
            'H5S_ALL',... in memory
            'H5S_ALL',... on disk
            'H5P_DEFAULT',... COLLECTIVE write change here for parallele
             field)  
  
  H5A.close ( attr_id )         
  H5S.close ( dataspace_id ) 
  H5D.close ( dataset_id )
  H5F.close ( file_id )
end