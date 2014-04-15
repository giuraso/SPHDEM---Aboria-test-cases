/*
 * species.h
 *
 * Copyright 2012 Martin Robinson
 *
  * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "BucketSort.h"

#include <vector>
//#include <boost/array.hpp>
//#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include "Zip.h"
#include "Vector.h"
//#include "MyRandom.h"

#ifndef HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#endif

namespace Aboria {

//#define DATA_typename ParticleInfo
//#define DATA_names   (r)(r0)(alive)(id)(saved_index)
//#define DATA_types   (Vect3d)(Vect3d)(bool)(int)(int)
//#include "Data.h"

template<typename DataType>
class Particles {
	template<typename T>
	friend class Particles;
public:
	class Value {
		friend class Particles;
	public:
		Value(){}
		Value(const Value& rhs) {
			r = rhs.r;
			r0 = rhs.r0;
			alive = rhs.alive;
			id = rhs.id;
			saved_index = rhs.saved_index;
			data = rhs.data;

			std::cout <<"copy constructor!!"<<std::endl;
		}
		Value& operator=(const Value &rhs) {
			if (this != &rhs) {
				r = rhs.r;
				r0 = rhs.r0;
				alive = rhs.alive;
				id = rhs.id;
				saved_index = rhs.saved_index;
				data = rhs.data;
			}
			std::cout <<"copying!!"<<std::endl;
			return *this;
		}
		const Vect3d& get_position() const {
			return r;
		}
		const Vect3d& get_old_position() const {
			return r0;
		}
		const DataType& get_data() const {
			return data;
		}
		DataType& get_data() {
			return data;
		}
		bool is_alive() const {
			return alive;
		}
		const std::size_t& get_id() const {
			return id;
		}
		const std::size_t& get_saved_index() const {
			return saved_index;
		}
		void mark_for_deletion() {
			alive = false;
		}
		template<typename T>
		boost::iterator_range<typename T::element_type::NeighbourSearch_type::const_iterator> get_neighbours(const T particles) {
			return boost::make_iterator_range(
			 particles->neighbour_search.find_broadphase_neighbours(get_position(), index,false),
			 particles->neighbour_search.end());
		}
		template<typename T>
		Vect3d correct_position_for_periodicity(const T particles, const Vect3d& position) {
			return particles->neighbour_search.correct_position_for_periodicity(r, position);
		}
	private:
		Vect3d r,r0;
		bool alive;
		std::size_t id;
		std::size_t index,saved_index;
		std::vector<size_t> neighbour_indicies;
		DataType data;
	};

	typedef typename std::vector<Value> data_type;
	typedef typename data_type::iterator iterator;
	typedef typename data_type::const_iterator const_iterator;
	struct get_pos {
		const Vect3d& operator()(const Value& i) const {
			return i.get_position();
		}
	};
	typedef BucketSort<const_iterator,get_pos> NeighbourSearch_type;


	Particles():
		next_id(0),
		neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
		searchable(false)
	{}

	static ptr<Particles<DataType> > New() {
		return ptr<Particles<DataType> >(new Particles<DataType> ());
	}

	Value& operator[](std::size_t idx) {
		return data[idx];
	}
	const Value& operator[](std::size_t idx) const {
		return const_cast<Value>(*this)[idx];
	}
	iterator begin() {
		return data.begin();
	}
	iterator end() {
		return data.end();
	}

//	void fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N) {
//		//TODO: assumes a 3d rectangular region
//		boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, boost::uniform_real<>(0,1));
//		const Vect3d dist = high-low;
//		for(int i=0;i<N;i++) {
//			add_particle(Vect3d(uni()*dist[0],uni()*dist[1],uni()*dist[2])+low);
//		}
//	}

	void delete_particles() {
		data.erase (std::remove_if( std::begin(data), std::end(data),
				[](Value& p) { return !p.is_alive(); }),
				std::end(data)
		);
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());
	}
	void clear() {
		data.clear();
	}
	size_t size() const {
		return data.size();
	}

	void save_indicies() {
		const int n = data.size();
		for (int i = 0; i < n; ++i) {
			data[i].saved_index = i;
		}
	}

	void init_neighbour_search(const Vect3d& low, const Vect3d& high, const double length_scale, const Vect3b& periodic) {
		neighbour_search.reset(low,high,length_scale,periodic);
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}
	

	template<typename F>
	void reset_neighbour_search(const double length_scale, F f) {
		std::for_each(begin(),end(),[&f](Value& i) {
			i.r0 = i.r;
			i.r = f(i);
		});
		neighbour_search.reset(neighbour_search.get_low(),neighbour_search.get_high(),length_scale,neighbour_search.get_periodic());
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}

	void reset_neighbour_search(const double length_scale) {
		neighbour_search.reset(neighbour_search.get_low(),neighbour_search.get_high(),length_scale,neighbour_search.get_periodic());
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}


	template<typename F>
	void create_particles(const int n, F f) {
		data.resize(data.size()+n);
		std::for_each(end()-n,end(),[this,&f](Value& i) {
			i.id = this->next_id++;
			i.alive = true;
			i.r = f(i);
			i.r0 = i.r;
		});
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());

	}

	template<typename F>
	void create_particles_cylinder(const Vect3d& min, const Vect3d& max, const int n, double dem_diameter, F f) {
		int shift=0;
		bool shiftz=false;
		double radius=(max[0]-min[0])/2-dem_diameter/2;
		const Vect3d origin(radius,radius,dem_diameter/2);
		int jmax=radius/(dem_diameter); //number of particles on radius
		for (int i = 0; i < n; ++i) { //z index
			radius=(max[0]-min[0])/2-dem_diameter/2;
			for (int j = 0; j < jmax; ++j) { //radius index
				radius=radius-dem_diameter;
				int kmax= radius*2*PI/dem_diameter; //number of particles on circumference
				
				double angle=shift*30+180*shiftz; // shift and shiftz avoid that vacuum is concentrated in certain point

						int index = data.size();
						data.resize(data.size()+kmax);
				for (int k = 0; k < kmax; ++k) {
					double d_angle=dem_diameter/radius;
					angle=angle+d_angle;
					data[index].r[0] = origin[0] + radius*cos(angle);
					data[index].r[1] = origin[1] + radius*sin(angle);
					data[index].r[2] = origin[2] + dem_diameter*i;
					data[index].r0 = data[index].r;
					data[index].id = next_id++;
					data[index].alive = true;
					f(data[index]);
					index++;
				}
				shift++;
			}
			shiftz=!shiftz;
		}
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());
	}


	template<typename F>
	void create_particles_grid(const Vect3d& startpoint, const Vect3d& endpoint, const Vect3i& n, F f) {
		const int nparticles = n.prod();
		int index = data.size();
		data.resize(data.size()+nparticles);
		Vect3d dx((endpoint[0]-startpoint[0])/n[0],(endpoint[1]-startpoint[1])/n[1],(endpoint[2]-startpoint[2])/n[2]);
		Vect3d low(startpoint[0]+dx[0]/2,startpoint[1]+dx[1]/2,startpoint[2]+dx[2]/2);
		for (int i = 0; i < n[0]; ++i) {
			for (int j = 0; j < n[1]; ++j) {
				for (int k = 0; k < n[2]; ++k) {

					data[index].r[0] = low[0] + dx[0]*(i);
					data[index].r[1] = low[1] + dx[1]*(j);
					data[index].r[2] = low[2] + dx[2]*(k);
					data[index].r0 = data[index].r;

					data[index].id = next_id++;
					data[index].alive = true;
					f(data[index]);
					index++;
				}
			}
		}
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());

	}

	template<typename F>
	void update_positions(iterator b, iterator e, F f) {
		std::for_each(b,e,[&f](Value& i) {
			i.r0 = i.r;
			i.r = f(i);
		});
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}
	template<typename F>
	void update_positions(F f) {
		std::for_each(begin(),end(),[&f](Value& i) {
			i.r0 = i.r;
			i.r = f(i);
		});
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}

#ifndef HAVE_VTK
	struct save_elem {
		typedef int result_type;
		save_elem(int index, vtkSmartPointer<vtkFloatArray>* datas):
			index(index),datas(datas){}
		template<typename T>
		result_type operator()(result_type state, T& t) const {
			datas[state]->SetValue(index,t);
			return state + 1;
		}
		//template<>
		result_type operator()(result_type state, Vect3d& t) const {
			datas[state]->SetTuple3(index,t[0],t[1],t[2]);
			return state + 1;
		}
		result_type operator()(result_type state, const Vect3d& t) const {
			datas[state]->SetTuple3(index,t[0],t[1],t[2]);
			return state + 1;
		}
		int index;
		vtkSmartPointer<vtkFloatArray>* datas;
	};
	struct setup_data {
		typedef int result_type;
		setup_data(vtkSmartPointer<vtkFloatArray>* datas):
			datas(datas){}

		template<typename T>
		result_type operator()(result_type state, T& t) const {
			datas[state]->SetNumberOfComponents(1);
			return state + 1;
		}
		result_type operator()(result_type state, const Vect3d& t) const {
			datas[state]->SetNumberOfComponents(3);
			return state + 1;
		}
		result_type operator()(result_type state, Vect3d& t) const {
			datas[state]->SetNumberOfComponents(3);
			return state + 1;
		}
		//template<>

		vtkSmartPointer<vtkFloatArray>* datas;
	};

	void  copy_to_vtk_grid(vtkSmartPointer<vtkUnstructuredGrid> grid) {
		vtkSmartPointer<vtkPoints> points = grid->GetPoints();
		if (!points) {
			points = vtkSmartPointer<vtkPoints>::New();
			grid->SetPoints(points);
		}
		vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
		if (!cells) {
			cells = vtkSmartPointer<vtkCellArray>::New();
			grid->SetCells(1,cells);
		}
		vtkSmartPointer<vtkUnsignedCharArray> cell_types = grid->GetCellTypesArray();
		vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(grid->GetPointData()->GetArray("id"));
		if (!ids) {
			ids = vtkSmartPointer<vtkIntArray>::New();
			ids->SetName("id");
			grid->GetPointData()->AddArray(ids);
		}
		constexpr size_t dn = std::tuple_size<DataType>::value;
		vtkSmartPointer<vtkFloatArray> datas[dn];
		for (int i = 0; i < dn; ++i) {
			char buffer[10];
			sprintf(buffer,"data%3d",i);
			datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(buffer));
			if (!datas[i]) {
				datas[i] = vtkSmartPointer<vtkFloatArray>::New();
				datas[i]->SetName(buffer);
				grid->GetPointData()->AddArray(datas[i]);
			}
		}
		boost::fusion::fold(DataType(),0,setup_data(datas));

		const vtkIdType n = size();
		points->SetNumberOfPoints(n);
		cells->Reset();
		cell_types->Reset();
		ids->SetNumberOfTuples(n);
		for (int i = 0; i < dn; ++i) {
			datas[i]->SetNumberOfTuples(n);
		}
		int j = 0;

		for(auto& i: *this) {
			const int index = j++;
			//std::cout <<"copying point at "<<i.get_position()<<" with id = "<<i.get_id()<<std::endl;
			points->SetPoint(index,i.get_position()[0],i.get_position()[1],i.get_position()[2]);
			cells->InsertNextCell(1);
			cells->InsertCellPoint(index);
			cell_types->InsertNextTuple1(1);
			ids->SetValue(index,i.get_id());


			boost::fusion::fold(i.get_data(),0,save_elem(index,datas));

//			auto v = std::get<0>(i.get_data());datas[0]->SetTuple3(index,v[0],v[1],v[2]);
//			v = std::get<1>(i.get_data());datas[1]->SetTuple3(index,v[0],v[1],v[2]);
//			v = std::get<2>(i.get_data());datas[2]->SetTuple3(index,v[0],v[1],v[2]);
//			v = std::get<3>(i.get_data());datas[3]->SetTuple3(index,v[0],v[1],v[2]);
//			double f = std::get<4>(i.get_data());datas[4]->SetValue(index,f);
		}

		//points->ComputeBounds();
	}
#endif

private:

	void enforce_domain(const Vect3d& low, const Vect3d& high, const Vect3b& periodic) {
		std::for_each(begin(),end(),[low,high,periodic](Value& i) {
			for (int d = 0; d < 3; ++d) {
				if (periodic[d]) {
					while (i.r[d]<low[d]) {
						i.r[d] += (high[d]-low[d]);
					}
					while (i.r[d]>=high[d]) {
						i.r[d] -= (high[d]-low[d]);
					}
				} else {
					if ((i.r[d]<low[d]) || (i.r[d]>=high[d])) {
						i.mark_for_deletion();
					}
				}
			}
		});
		if ((periodic[0]==false)||(periodic[1]==false)||(periodic[2]==false)) {
			data.erase (std::remove_if( std::begin(data), std::end(data),
					[](Value& p) { return !p.is_alive(); }),
					std::end(data)
			);
		}
	}


	data_type data;
	bool searchable;
	int next_id;
	NeighbourSearch_type neighbour_search;

};


}


#endif /* SPECIES_H_ */
