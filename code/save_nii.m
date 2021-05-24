%  Save NIFTI dataset. Support both *.nii and *.hdr/*.img file extension.
%  If file extension is not provided, *.hdr/*.img will be used as default.
%  
%  Usage: save_nii(nii, filename, [old_RGB])
%  
%  nii.hdr - struct with NIFTI header fields (from load_nii.m or make_nii.m)
%
%  nii.img - 3D (or 4D) matrix of NIFTI data.
%
%  filename - NIFTI file name.
%
%  old_RGB    - an optional boolean variable to handle special RGB data 
%       sequence [R1 R2 ... G1 G2 ... B1 B2 ...] that is used only by 
%       AnalyzeDirect (Analyze Software). Since both NIfTI and Analyze
%       file format use RGB triple [R1 G1 B1 R2 G2 B2 ...] sequentially
%       for each voxel, this variable is set to FALSE by default. If you
%       would like the saved image only to be opened by AnalyzeDirect 
%       Software, set old_RGB to TRUE (or 1). It will be set to 0, if it
%       is default or empty.
%  
%  Tip: to change the data type, set nii.hdr.dime.datatype,
%	and nii.hdr.dime.bitpix to:
% 
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
%   128 uint RGB                  (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
%   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  
%  Part of this file is copied and modified from:
%  http://www.mathworks.com/matlabcentral/fileexchange/1878-mri-analyze-tools
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%  - "old_RGB" related codes in "save_nii.m" are added by Mike Harms (2006.06.28) 
%
function save_nii(nii, fileprefix, old_RGB)
   
   if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
	~isfield(nii,'img') | ~exist('fileprefix','var') | isempty(fileprefix)

      error('Usage: save_nii(nii, filename, [old_RGB])');
   end

   if isfield(nii,'untouch') & nii.untouch == 1
      error('Usage: please use ''save_untouch_nii.m'' for the untouched structure.');
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end
   
   filetype = 1;

   %  Note: fileprefix is actually the filename you want to save
   %   
   if findstr('.nii',fileprefix)
      filetype = 2;
      fileprefix = strrep(fileprefix,'.nii','');
   end
   
   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end
   
   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end
   
   write_nii(nii, filetype, fileprefix, old_RGB);

   if filetype == 1

      %  So earlier versions of SPM can also open it with correct originator
      %
      M=[[diag(nii.hdr.dime.pixdim(2:4)) -[nii.hdr.hist.originator(1:3).*nii.hdr.dime.pixdim(2:4)]'];[0 0 0 1]];
      save([fileprefix '.mat'], 'M');
   end
   
   return					% save_nii


%-----------------------------------------------------------------------------------
function write_nii(nii, filetype, fileprefix, old_RGB)

   hdr = nii.hdr;

   if isfield(nii,'ext') & ~isempty(nii.ext)
      ext = nii.ext;
      [ext, esize_total] = verify_nii_ext(ext);
   else
      ext = [];
   end

   switch double(hdr.dime.datatype),
   case   1,
      hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = int16(8 ); precision = 'uint8';
   case   4,
      hdr.dime.bitpix = int16(16); precision = 'int16';
   case   8,
      hdr.dime.bitpix = int16(32); precision = 'int32';
   case  16,
      hdr.dime.bitpix = int16(32); precision = 'float32';
   case  32,
      hdr.dime.bitpix = int16(64); precision = 'float32';
   case  64,
      hdr.dime.bitpix = int16(64); precision = 'float64';
   case 128,
      hdr.dime.bitpix = int16(24); precision = 'uint8';
   case 256 
      hdr.dime.bitpix = int16(8 ); precision = 'int8';
   case 511,
      hdr.dime.bitpix = int16(96); precision = 'float32';
   case 512 
      hdr.dime.bitpix = int16(16); precision = 'uint16';
   case 768 
      hdr.dime.bitpix = int16(32); precision = 'uint32';
   case 1024
      hdr.dime.bitpix = int16(64); precision = 'int64';
   case 1280
      hdr.dime.bitpix = int16(64); precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = int16(128); precision = 'float64';
   otherwise
      error('This datatype is not supported');
   end
   
   hdr.dime.glmax = round(double(max(nii.img(:))));
   hdr.dime.glmin = round(double(min(nii.img(:))));
   
   if filetype == 2
      fid = fopen(sprintf('%s.nii',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.nii.',fileprefix);
         error(msg);
      end
      
      hdr.dime.vox_offset = 352;

      if ~isempty(ext)
         hdr.dime.vox_offset = hdr.dime.vox_offset + esize_total;
      end

      hdr.hist.magic = 'n+1';
      save_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end
   else
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');
      
      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end
      
      hdr.dime.vox_offset = 0;
      hdr.hist.magic = 'ni1';
      save_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end
      
      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   end

   ScanDim = double(hdr.dime.dim(5));		% t
   SliceDim = double(hdr.dime.dim(4));		% z
   RowDim   = double(hdr.dime.dim(3));		% y
   PixelDim = double(hdr.dime.dim(2));		% x
   SliceSz  = double(hdr.dime.pixdim(4));
   RowSz    = double(hdr.dime.pixdim(3));
   PixelSz  = double(hdr.dime.pixdim(2));
   
   x = 1:PixelDim;

   if filetype == 2 & isempty(ext)
      skip_bytes = double(hdr.dime.vox_offset) - 348;
   else
      skip_bytes = 0;
   end

   if double(hdr.dime.datatype) == 128

      %  RGB planes are expected to be in the 4th dimension of nii.img
      %
      if(size(nii.img,4)~=3)
         error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
      end

      if old_RGB
         nii.img = permute(nii.img, [1 2 4 3 5 6 7 8]);
      else
         nii.img = permute(nii.img, [4 1 2 3 5 6 7 8]);
      end
   end

   if double(hdr.dime.datatype) == 511

      %  RGB planes are expected to be in the 4th dimension of nii.img
      %
      if(size(nii.img,4)~=3)
         error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
      end

      if old_RGB
         nii.img = permute(nii.img, [1 2 4 3 5 6 7 8]);
      else
         nii.img = permute(nii.img, [4 1 2 3 5 6 7 8]);
      end
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      real_img = real(nii.img(:))';
      nii.img = imag(nii.img(:))';
      nii.img = [real_img; nii.img];
   end

   if skip_bytes
      fwrite(fid, zeros(1,skip_bytes), 'uint8');
   end

   fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
   fclose(fid);

   return;					% write_nii

%  Save NIFTI header extension.
%
%  Usage: save_nii_ext(ext, fid)
%
%  ext - struct with NIFTI header extension fields.
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function save_nii_ext(ext, fid)

   if ~exist('ext','var') | ~exist('fid','var')
      error('Usage: save_nii_ext(ext, fid)');
   end

   if ~isfield(ext,'extension') | ~isfield(ext,'section') | ~isfield(ext,'num_ext')
      error('Wrong header extension');
   end

   write_ext(ext, fid);

   return;                                      % save_nii_ext


%---------------------------------------------------------------------
function write_ext(ext, fid)

   fwrite(fid, ext.extension, 'uchar');

   for i=1:ext.num_ext
      fwrite(fid, ext.section(i).esize, 'int32');
      fwrite(fid, ext.section(i).ecode, 'int32');
      fwrite(fid, ext.section(i).edata, 'uchar');
   end

   return;                                      % write_ext

%  internal function
  
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function save_nii_hdr(hdr, fid)
   
   if ~exist('hdr','var') | ~exist('fid','var')
      error('Usage: save_nii_hdr(hdr, fid)');
   end
   
   if ~isequal(hdr.hk.sizeof_hdr,348),
      error('hdr.hk.sizeof_hdr must be 348.');
   end
   
   if hdr.hist.qform_code == 0 & hdr.hist.sform_code == 0
      hdr.hist.sform_code = 1;
      hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
      hdr.hist.srow_x(2) = 0;
      hdr.hist.srow_x(3) = 0;
      hdr.hist.srow_y(1) = 0;
      hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
      hdr.hist.srow_y(3) = 0;
      hdr.hist.srow_z(1) = 0;
      hdr.hist.srow_z(2) = 0;
      hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
      hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
      hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
      hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
   end
   
   write_header(hdr, fid);

   return;					% save_nii_hdr


%---------------------------------------------------------------------
function write_header(hdr, fid)

        %  Original header structures
	%  struct dsr				/* dsr = hdr */
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/
   
   header_key(fid, hdr.hk);
   image_dimension(fid, hdr.dime);
   data_history(fid, hdr.hist);
   
   %  check the file size is 348 bytes
   %
   fbytes = ftell(fid);
   
   if ~isequal(fbytes,348),
      msg = sprintf('Header size is not 348 bytes.');
      warning(msg);
   end
    
   return;					% write_header


%---------------------------------------------------------------------
function header_key(fid, hk)
   
   fseek(fid,0,'bof');

	%  Original header structures    
	%  struct header_key                      /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
        
   fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.
    
   % data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
   % fwrite(fid, data_type(1:10), 'uchar');
   pad = zeros(1, 10-length(hk.data_type));
   hk.data_type = [hk.data_type  char(pad)];
   fwrite(fid, hk.data_type(1:10), 'uchar');
    
   % db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
   % fwrite(fid, db_name(1:18), 'uchar');
   pad = zeros(1, 18-length(hk.db_name));
   hk.db_name = [hk.db_name  char(pad)];
   fwrite(fid, hk.db_name(1:18), 'uchar');
    
   fwrite(fid, hk.extents(1),       'int32');
   fwrite(fid, hk.session_error(1), 'int16');
   fwrite(fid, hk.regular(1),       'uchar');	% might be uint8
    
   % fwrite(fid, hk.hkey_un0(1),    'uchar');
   % fwrite(fid, hk.hkey_un0(1),    'uint8');
   fwrite(fid, hk.dim_info(1),      'uchar');
    
   return;					% header_key


%---------------------------------------------------------------------
function image_dimension(fid, dime)

	%  Original header structures        
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%				pixdim[4] - volume timing, in msec
	%					..etc
	%			*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
   fwrite(fid, dime.dim(1:8),        'int16');
   fwrite(fid, dime.intent_p1(1),  'float32');
   fwrite(fid, dime.intent_p2(1),  'float32');
   fwrite(fid, dime.intent_p3(1),  'float32');
   fwrite(fid, dime.intent_code(1),  'int16');
   fwrite(fid, dime.datatype(1),     'int16');
   fwrite(fid, dime.bitpix(1),       'int16');
   fwrite(fid, dime.slice_start(1),  'int16');
   fwrite(fid, dime.pixdim(1:8),   'float32');
   fwrite(fid, dime.vox_offset(1), 'float32');
   fwrite(fid, dime.scl_slope(1),  'float32');
   fwrite(fid, dime.scl_inter(1),  'float32');
   fwrite(fid, dime.slice_end(1),    'int16');
   fwrite(fid, dime.slice_code(1),   'uchar');
   fwrite(fid, dime.xyzt_units(1),   'uchar');
   fwrite(fid, dime.cal_max(1),    'float32');
   fwrite(fid, dime.cal_min(1),    'float32');
   fwrite(fid, dime.slice_duration(1), 'float32');
   fwrite(fid, dime.toffset(1),    'float32');
   fwrite(fid, dime.glmax(1),        'int32');
   fwrite(fid, dime.glmin(1),        'int32');
   
   return;					% image_dimension


%---------------------------------------------------------------------
function data_history(fid, hist)
    
	% Original header structures
	%struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */
	
   % descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
   % fwrite(fid, descrip(1:80),    'uchar');
   pad = zeros(1, 80-length(hist.descrip));
   hist.descrip = [hist.descrip  char(pad)];
   fwrite(fid, hist.descrip(1:80), 'uchar');
    
   % aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
   % fwrite(fid, aux_file(1:24),   'uchar');
   pad = zeros(1, 24-length(hist.aux_file));
   hist.aux_file = [hist.aux_file  char(pad)];
   fwrite(fid, hist.aux_file(1:24), 'uchar');
    
   fwrite(fid, hist.qform_code,    'int16');
   fwrite(fid, hist.sform_code,    'int16');
   fwrite(fid, hist.quatern_b,   'float32');
   fwrite(fid, hist.quatern_c,   'float32');
   fwrite(fid, hist.quatern_d,   'float32');
   fwrite(fid, hist.qoffset_x,   'float32');
   fwrite(fid, hist.qoffset_y,   'float32');
   fwrite(fid, hist.qoffset_z,   'float32');
   fwrite(fid, hist.srow_x(1:4), 'float32');
   fwrite(fid, hist.srow_y(1:4), 'float32');
   fwrite(fid, hist.srow_z(1:4), 'float32');

   % intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
   % fwrite(fid, intent_name(1:16),    'uchar');
   pad = zeros(1, 16-length(hist.intent_name));
   hist.intent_name = [hist.intent_name  char(pad)];
   fwrite(fid, hist.intent_name(1:16), 'uchar');
    
   % magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
   % fwrite(fid, magic(1:4),           'uchar');
   pad = zeros(1, 4-length(hist.magic));
   hist.magic = [hist.magic  char(pad)];
   fwrite(fid, hist.magic(1:4),        'uchar');
    
   return;					% data_history

