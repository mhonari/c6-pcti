%%This code builds a 6 pillar model and runs a non-parametric study and export the eigenfrequencies to file. 
%%created with COMSOL 5.3

import com.comsol.model.*
import com.comsol.model.util.*

LATTICE_CONST = 300;

model = ModelUtil.create('Model');

model.modelPath('F:\Mostafa\comsol\clean model\export_test');

model.component.create('comp1', true);



%%parameters
model.param.set('k', '1', 'fraction of k amplitude. Varies between 0 and 1');
model.param.set('a1y', '0[m]', 'First lattice vector, y component');
model.param.set('eps_a', '1', 'dielectric const of the surrounding');
model.param.set('band', '1', 'band number');
model.param.set('c', '1', 'normalized speed of light');
model.param.set('eps_d', '12', 'dielectric constant of pillars');
model.param.set('c1', '2.9', 'geometry param1');
model.param.set('c2', '3', 'geometry param2');
model.param.set('a', strcat(num2str(LATTICE_CONST, 15), '[nm]') , 'Lattice constant');
model.param.set('a1x', 'a', 'First lattice vector, x component');
model.param.set('a2x', 'a/2', 'Second lattice vector, x component');
model.param.set('R', 'a/c1', 'Pillar center distance');
model.param.set('L', 'a/sqrt(3)', 'Hexagon side');
model.param.set('a2y', 'a*sqrt(3)/2', 'Second lattice vector, y component');
model.param.set('d', 'R/c2', 'pillars'' radii');

model.component('comp1').variable.create('var1');

model.component('comp1').variable('var1').set('kx_KM', '(-(k-2)+4)*pi/(3*a)');
model.component('comp1').variable('var1').descr('kx_KM', 'kx for KM path (2<k<3)');
model.component('comp1').variable('var1').set('kx_MG', '(1-k)*0.5*(b1x)');
model.component('comp1').variable('var1').descr('kx_MG', 'kx for M-Gamma path');
model.component('comp1').variable('var1').set('ky_MG', '-(1-k)*0.5*(b1y)');
model.component('comp1').variable('var1').descr('ky_MG', 'ky for M-Gamma path');
model.component('comp1').variable('var1').set('ky_KM', '(k-2)*pi/(sqrt(3)*a)');
model.component('comp1').variable('var1').descr('ky_KM', 'ky for KM path (2<k<3)');
model.component('comp1').variable('var1').set('ky_GK', '0[1/m]');
model.component('comp1').variable('var1').descr('ky_GK', 'ky for Gamma-K path (1<k<2)');
model.component('comp1').variable('var1').set('b2y', '2*pi*a1x/(a1x*a2y-a1y*a2x)');
model.component('comp1').variable('var1').descr('b2y', 'Second reciprocal lattice vector, y component');
model.component('comp1').variable('var1').set('b2x', '-2*pi*a1y/(a1x*a2y-a1y*a2x)');
model.component('comp1').variable('var1').descr('b2x', 'Second reciprocal lattice vector, x component');
model.component('comp1').variable('var1').set('b1y', '-2*pi*a2x/(a1x*a2y-a1y*a2x)');
model.component('comp1').variable('var1').descr('b1y', 'First reciprocal lattice vector, y component');
model.component('comp1').variable('var1').set('b1x', '2*pi*a2y/(a1x*a2y-a1y*a2x)');
model.component('comp1').variable('var1').descr('b1x', 'First reciprocal lattice vector, x component');
model.component('comp1').variable('var1').set('kx_GK', '4*pi/(3*a)*(k-1)');
model.component('comp1').variable('var1').descr('kx_GK', 'kx for Gamma-K path (1<k<2)');
model.component('comp1').variable('var1').set('kxx', 'if( k<1 , kx_MG  ,  if( k<2  ,  kx_GK   , kx_KM)  )');
model.component('comp1').variable('var1').descr('kxx', 'kx for MGKM path');
model.component('comp1').variable('var1').set('kyy', 'if( k<1 , ky_MG  ,  if( k<2  ,  ky_GK   , ky_KM)  )');
model.component('comp1').variable('var1').descr('kyy', 'ky for MGKM path');

%%Geometry
model.component('comp1').geom.create('geom1', 2);
%hexagon
model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 0, 0, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 'L', 0, 1);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', '-L*sqrt(3)/2', 1, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 'L/2', 1, 1);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', '-L*sqrt(3)/2', 2, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', '-L/2', 2, 1);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 0, 3, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', '-L', 3, 1);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 'L*sqrt(3)/2', 4, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', '-L/2', 4, 1);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 'L*sqrt(3)/2', 5, 0);
model.component('comp1').geom('geom1').feature('pol1').setIndex('table', 'L/2', 5, 1);
model.component('comp1').geom('geom1').run('pol1');
model.geom('geom1').run('pol1');
%pillars
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('r', 'd');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'R' '0'});
model.component('comp1').geom('geom1').run('c1');
model.component('comp1').geom('geom1').create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('c2').set('r', 'd');
model.component('comp1').geom('geom1').feature('c2').set('pos', {'0' 'R*sqrt(3)/2'});
model.component('comp1').geom('geom1').feature('c2').setIndex('pos', 'R/2', 0);
model.component('comp1').geom('geom1').run('c2');
model.component('comp1').geom('geom1').create('c3', 'Circle');
model.component('comp1').geom('geom1').feature('c3').set('r', 'd');
model.component('comp1').geom('geom1').feature('c3').set('pos', {'0' 'R*sqrt(3)/2'});
model.component('comp1').geom('geom1').feature('c3').setIndex('pos', '-R/2', 0);
model.component('comp1').geom('geom1').run('c3');
model.component('comp1').geom('geom1').create('c4', 'Circle');
model.component('comp1').geom('geom1').feature('c4').set('r', 'd');
model.component('comp1').geom('geom1').feature('c4').set('pos', {'-R' '0'});
model.component('comp1').geom('geom1').run('c4');
model.component('comp1').geom('geom1').create('c5', 'Circle');
model.component('comp1').geom('geom1').feature('c5').set('pos', {'-R/2' '-R*sqrt(3)/2'});
model.component('comp1').geom('geom1').feature('c5').set('r', 'd');
model.component('comp1').geom('geom1').run('c5');
model.component('comp1').geom('geom1').create('c6', 'Circle');
model.component('comp1').geom('geom1').feature('c6').set('pos', {'R/2' '-R*sqrt(3)/2'});
model.component('comp1').geom('geom1').feature('c6').set('r', 'd');
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').run;

%%materials
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').selection.set([2 3 4 5 6 7]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'sqrt(eps_d)'});
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').selection.set([1]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'sqrt(eps_a)'});
model.component('comp1').material('mat1').label('rods');
model.component('comp1').material('mat2').label('background');

%%Physics
model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').prop('components').set('components', 'outofplane');
model.component('comp1').physics('ewfd').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('ewfd').feature('pc1').selection.set([1 6]);
model.component('comp1').physics('ewfd').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc1').set('kFloquet', {'kxx' 'kyy' '0'});
model.component('comp1').physics('ewfd').feature('pc1').set('constraintMethod', 'nodal');
model.component('comp1').physics('ewfd').create('pc2', 'PeriodicCondition', 1);
model.component('comp1').physics('ewfd').feature('pc2').selection.set([3 4]);
model.component('comp1').physics('ewfd').feature('pc2').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc2').set('kFloquet', {'kxx' 'kyy' '0'});
model.component('comp1').physics('ewfd').feature('pc2').set('constraintMethod', 'nodal');
model.component('comp1').physics('ewfd').create('pc3', 'PeriodicCondition', 1);
model.component('comp1').physics('ewfd').feature('pc3').selection.set([2 5]);
model.component('comp1').physics('ewfd').feature('pc3').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc3').set('kFloquet', {'kxx' 'kyy' '0'});
model.component('comp1').physics('ewfd').feature('pc3').set('constraintMethod', 'nodal');

%%Mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 2);
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([1 2 3]);
%copy edge
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([2]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([5]);
model.component('comp1').mesh('mesh1').create('cpe2', 'CopyEdge');
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').set([1]);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').set([6]);
model.component('comp1').mesh('mesh1').create('cpe3', 'CopyEdge');
model.component('comp1').mesh('mesh1').feature('cpe3').selection('source').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('destination').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('source').set([3]);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('destination').set([4]);
%free triangular mesh
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').run;

%%Study settings
model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').activate('ewfd', true);
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('neigs', 7);
model.study('std1').feature('eig').set('shiftactive', true);
model.study('std1').feature('eig').set('shift', 'c_const/a*0.1');


model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
model.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;


% model.result.create('pg2', 'PlotGroup2D');
% model.result('pg2').label('Electric Field (ewfd)');
% model.result('pg2').set('frametype', 'spatial');
% model.result('pg2').set('data', 'dset1');
% model.result('pg2').feature.create('surf1', 'Surface');
% model.result('pg2').feature('surf1').set('data', 'parent');
% model.result('pg2').run;

%%create 1d band diagram plot
model.result.create('pg3', 'PlotGroup1D');
model.result('pg3').create('glob1', 'Global');
model.result('pg3').feature('glob1').set('expr', {'ewfd.freq'});
model.result('pg3').feature('glob1').set('descr', {'Frequency'});
model.result('pg3').feature('glob1').set('unit', {'Hz'});
model.result('pg3').set('data', 'dset1');
model.result('pg3').feature('glob1').set('xdatasolnumtype', 'outer');
model.result('pg3').feature('glob1').set('xdata', 'expr');
model.result('pg3').feature('glob1').set('xdataexpr', 'k');
model.result('pg3').run;

%%Get the eigenfrequencies
sol_info = mphsolinfo(model);
sol_vals = sol_info.solvals;
eigfreqs = -imag(sol_vals)*300e-9/(2*pi*3e8);
% %export the eigenfrequencies
% model.result.export.create('plot1', 'pg3', 'glob1', 'Plot');
% model.result.export('plot1').set('filename', 'F:\Mostafa\comsol\clean model\export_test\1d_export__test.txt');
% model.result.export('plot1').set('header', false);
% model.result.export('plot1').run;