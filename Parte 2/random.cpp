/*
 * random.cpp
 * 
 * Copyright 2020 otavio <otavio@otavio-Aspire-5542>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <random>

int main(){
	
	std::uniform_real_distribution<double> unif(-1,1);
			std::default_random_engine re;
		for(int i=0;i<100;i++){
			double ruido = unif(re);
			std::cout<<ruido<<std::endl;		
		}
	return 0;
}

