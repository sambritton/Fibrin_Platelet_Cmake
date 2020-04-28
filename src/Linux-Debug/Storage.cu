#include <sys/stat.h>

#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "System.h"
#include "System_Builder.h"
#include "SystemStructures.h"
#include "Storage.h"


Storage::Storage(std::weak_ptr<System> a_system,
	std::weak_ptr<SystemBuilder> b_system ) {
	//std::cout << "FDM constructor" << std::endl;

	system = a_system;
	builder = b_system;

	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		std::stringstream stream_min;
		std::stringstream stream_max;

		unsigned domain_size = ceil((sys->domainParams.maxX + 
			sys->domainParams.maxY + 
			sys->domainParams.maxZ) / 3.0);

		unsigned plt_count = sys->generalParams.maxPltCount;
		unsigned num_filo = sys->generalParams.plt_tndrl_intrct;
		double max_force = sys->generalParams.max_dynamic_plt_force;
		double min_force = sys->generalParams.pltForce;
		bool response_nonlin = sys->generalParams.use_nonlinear_dynamic_force;
		bool response = sys->generalParams.use_dynamic_plt_force;
		bool force_scale = sys->generalParams.distribute_plt_force;
		
		stream_min << std::fixed << std::setprecision(2) << min_force;
		std::string str_min_force = stream_min.str();
		
		stream_max << std::fixed << std::setprecision(2) << max_force;
		std::string str_max_force = stream_max.str();

		std::string str_force_scale = "_pltForceScale_";
		std::string str_nonlinear_response = "_nonlinDynamicResp_";
		std::string str_dynamic_response = "_dynamicResp_";
		std::string str_domain = "_domain_";
		std::string str_plt_count = "_plt_count_";
		std::string str_filo_count = "_filo_count_";
		std::string str_maxF = "_maxForce_";
		std::string str_minF = "_minForce_";

		std::string str_a = "Animation_";
		std::string str_p = "Params_";
		std::string str_s = "State_";

		str_animation = str_a +
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		std::cout<<"domain: " << domain_size<<
			" str_plt_count: " << plt_count<<
			" str_filo_count: " << num_filo<<
			" str_minF: " << min_force<< 
			" str_response: " << response << std::endl;


		const int dir_err_anim = mkdir(str_animation.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_anim)
		{
			printf("Error creating directory animation test!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_animation.c_str());
		}

		str_params = str_p + 
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		const int dir_err_params = mkdir(str_params.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_params)
		{
			printf("Error creating directory params!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_params.c_str());
		}

		str_state = str_s + 
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		const int dir_err_state = mkdir(str_state.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_state)
		{
			printf("Error creating directory state!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_state.c_str());
		}
	}
};

void Storage::save_current_state(){
	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = str_state+"/State_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		//unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		
		//place nodes
		ofs << std::setprecision(30) <<std::fixed<< "time " << sys->generalParams.currentTime<<std::endl;
		
		//thrust::copy(sys->nodeInfoVecs.nodeLocX.begin(),
		//	sys->nodeInfoVecs.nodeLocX.end(), hostLocX.begin());
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			
			double x = sys->nodeInfoVecs.nodeLocX[i];
			double y = sys->nodeInfoVecs.nodeLocY[i];
			double z = sys->nodeInfoVecs.nodeLocZ[i];
			ofs << std::setprecision(30) <<std::fixed<< "node " << x << " " << y << " " << z <<std::endl;
		}
		
		//place plts
		for (unsigned i = 0; i < sys->pltInfoVecs.pltLocX.size(); i++) {
			double x = sys->pltInfoVecs.pltLocX[i];
			double y = sys->pltInfoVecs.pltLocY[i];
			double z = sys->pltInfoVecs.pltLocZ[i];
			ofs << std::setprecision(30) <<std::fixed<< "plt " << x << " " << y << " " << z <<std::endl;
		
		}
		
		//
		ofs << std::setprecision(30) <<std::fixed<< "edge_count " << sys->generalParams.currentEdgeCount<<std::endl;
		for (unsigned i = 0; i < sys->generalParams.currentEdgeCount; i++) {
			unsigned left = sys->nodeInfoVecs.hostEdgeLeft[i];
			unsigned right = sys->nodeInfoVecs.hostEdgeRight[i];

			ofs << std::setprecision(30) <<std::fixed<< "edge_lr " << left << " " << right <<std::endl;
		}

		/*
		ofs << std::setprecision(30) <<std::fixed<< "minX " << sys->domainParams.originMinX<<std::endl;
		ofs << std::setprecision(30) <<std::fixed<< "maxX " << sys->domainParams.originMaxX<<std::endl;
		ofs << std::setprecision(30) <<std::fixed<< "minY " << sys->domainParams.originMinY<<std::endl;
		ofs << std::setprecision(30) <<std::fixed<< "maxY " << sys->domainParams.originMaxY<<std::endl;
		ofs << std::setprecision(30) <<std::fixed<< "minZ " << sys->domainParams.originMinZ<<std::endl;
		ofs << std::setprecision(30) <<std::fixed<< "maxZ " << sys->domainParams.originMaxZ<<std::endl;
		*/
		for (unsigned i = 0; i < sys->wlcInfoVecs.globalNeighbors.size(); i++) {
			unsigned nbr = sys->wlcInfoVecs.globalNeighbors[i];

			ofs << std::setprecision(30) <<std::fixed<< "global_nbr " << nbr << std::endl;
		}

		for (unsigned i = 0; i < sys->wlcInfoVecs.lengthZero.size(); i++) {
			double len_zero = sys->wlcInfoVecs.lengthZero[i];

			ofs << std::setprecision(30) <<std::fixed<< "length_zero " << len_zero << std::endl;
		}

		for (unsigned i = 0; i < sys->nodeInfoVecs.isNodeInPltVol.size(); i++) {
			bool is_in_plt = sys->nodeInfoVecs.isNodeInPltVol[i];

			ofs << std::setprecision(3) <<std::fixed<< "is_node_in_plt " << is_in_plt << std::endl;
		}
		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeId.size(); i++) {
			unsigned fil_node_id = sys->pltInfoVecs.tndrlNodeId[i];

			ofs << std::setprecision(3) <<std::fixed<< "file_node_id " << fil_node_id << std::endl;
		}
		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeType.size(); i++) {
			unsigned fil_node_type = sys->pltInfoVecs.tndrlNodeType[i];
			ofs << std::setprecision(3) <<std::fixed<< "fil_node_type " << fil_node_type << std::endl;
		}

		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeNum.size(); i++) {
			unsigned num_filopodia = sys->pltInfoVecs.tndrlNodeNum[i];
			ofs << std::setprecision(3) <<std::fixed<< "filopodia_count " << num_filopodia << std::endl;
		}
		
	}
}

void Storage::save_params(void) {
	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = str_params+"/Param_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		//unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "time " << sys->generalParams.currentTime<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minX " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxX " << sys->domainParams.maxX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minY " << sys->domainParams.minY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxY " << sys->domainParams.maxY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minZ " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxZ " << sys->domainParams.maxX<<std::endl;

		ofs << std::setprecision(5) <<std::fixed<< "total_applied_force_fibers " << sys->generalParams.totalAppliedForce_Fiber<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "total_applied_force_plts " << sys->generalParams.totalAppliedForce_Plt<<std::endl;

		ofs << std::setprecision(5) <<std::fixed<< "original_node_count " << originalNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "node_count_discretize " << maxNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_edge_count " << originalEdgeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "edge_count_discretize " << edgeCountDiscretize <<std::endl;
		
		//place nodes
		
		//thrust::copy(sys->nodeInfoVecs.nodeLocX.begin(),hostEdgeRight
		//	sys->nodeInfoVecs.nodeLocX.end(), hostLocX.begin());
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			
			double x = sys->nodeInfoVecs.nodeLocX[i];
			double y = sys->nodeInfoVecs.nodeLocY[i];
			double z = sys->nodeInfoVecs.nodeLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "node " << x << " " << y << " " << z <<std::endl;
		}
		
		//place plts
		for (unsigned i = 0; i < sys->pltInfoVecs.pltLocX.size(); i++) {
			double x = sys->pltInfoVecs.pltLocX[i];
			double y = sys->pltInfoVecs.pltLocY[i];
			double z = sys->pltInfoVecs.pltLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "plt " << x << " " << y << " " << z <<std::endl;
		
		}

				
		//place plts force
		for (unsigned i = 0; i < sys->pltInfoVecs.pltForceX.size(); i++) {
			double x = sys->pltInfoVecs.pltForceX[i];
			double y = sys->pltInfoVecs.pltForceY[i];
			double z = sys->pltInfoVecs.pltForceZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "force_on_plt " <<
			x << " " << y << " " << z <<std::endl;
		}

		//now print plt and corresponding node it's connected to		
		unsigned maxPltCount = sys->pltInfoVecs.pltLocX.size();
		unsigned num_connections = sys->pltInfoVecs.numConnections;
		for (unsigned edge = 0; edge < num_connections; edge++ ){

			//because nodes are placed after platelets, their vtk file id is incremented. 
			//notice that this represents the vtk id, not the id within the c++ program
			unsigned node_id_vtk = maxPltCount + edge;
			unsigned node_id = sys->pltInfoVecs.nodeImagingConnection[edge];
			
			unsigned plt_id = sys->pltInfoVecs.pltImagingConnection[edge];
			
			ofs << std::setprecision(5) << std::fixed<< 
			"plt_id_node_id "<< plt_id << " "<< node_id_vtk <<std::endl;

		}


		//place force node is experiencing (vector form)
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			double f_x = sys->nodeInfoVecs.nodeForceX[i];
			double f_y = sys->nodeInfoVecs.nodeForceY[i];
			double f_z = sys->nodeInfoVecs.nodeForceZ[i];
			
			ofs << std::setprecision(5) <<std::fixed<< "force_on_node " << 
			f_x <<" " << f_y <<" " << f_z << std::endl;
		
		}

		//place original edges
		for (unsigned edge = 0; edge < sys->generalParams.originEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];
			ofs <<"original_edge_discretized " <<idL <<" "<< idR <<std::endl;
			
		}
				 
		//place added edges
		for (unsigned edge = sys->generalParams.originEdgeCount; edge < sys->generalParams.currentEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];
			ofs <<"added_edge " <<idL <<" "<< idR <<std::endl;
			
		}

		//original edge strain
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];

			ofs << std::setprecision(5)<< std::fixed<<"original_edge_strain " << val <<std::endl;
		}
				
		//original edge alignment
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeAlignment[i];
			ofs << std::setprecision(5)<< std::fixed<<"original_edge_alignment " << val <<std::endl;
		}

		//added edge strain
		for (unsigned i = sys->generalParams.originEdgeCount; i < sys->generalParams.currentEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];
			ofs << std::setprecision(5)<< std::fixed<<"added_edge_strain " << val <<std::endl;
		}
		
		//added links per node.
		for (unsigned i = 0; i < sys->generalParams.maxNodeCount; i++ ){
			unsigned val = sys->wlcInfoVecs.currentNodeEdgeCountVector[i] - 
				sys->wlcInfoVecs.numOriginalNeighborsNodeVector[i];
			ofs << std::setprecision(5)<< std::fixed<<"bind_sites_per_node " << val <<std::endl;
		}



	}
};


void Storage::print_VTK_File() {

	std::shared_ptr<System> sys = system.lock();
	if (sys) {
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation + "/FibrinNetwork_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());


		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned maxNeighborCount = (sys->generalParams).maxNeighborCount;

		unsigned numEdges = sys->generalParams.currentEdgeCount;

		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ofs << "POINTS " << maxNodeCount << " float" << std::endl;
		for (unsigned i = 0; i< maxNodeCount; i++) { 
			double xPos = sys->nodeInfoVecs.nodeLocX[i];
			double yPos = sys->nodeInfoVecs.nodeLocY[i];
			double zPos = sys->nodeInfoVecs.nodeLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles


		unsigned numCells = numEdges;
		unsigned numNumsInCells = 3 * numEdges;


		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;

		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];

			ofs<< 2 << " " << idL << " " << idR << std::endl;
		}
	/*	for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->wlcInfoVecs.globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX

				//counts only half
				if ((idA < idB) && (idB < maxNodeCount) ) {
					ofs<< 2 << " " << idA << " " << idB << std::endl;
				}
			}
		}*/

		ofs << "CELL_TYPES " << numCells << std::endl;
		for (unsigned i = 0; i< numEdges; i++) {
			ofs << 3 << std::endl; //edge joining two points
		}



		//
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Fiber_Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idB = sys->nodeInfoVecs.hostEdgeRight[edge];

			unsigned begin = idA * sys->generalParams.maxNeighborCount;
			unsigned end = begin + sys->generalParams.maxNeighborCount;
			double L0;
			for (unsigned i = begin; i < end; i++) {
				unsigned idTemp = sys->wlcInfoVecs.globalNeighbors[i];
				if (idTemp == idB){
					L0 = sys->wlcInfoVecs.lengthZero[i];
				}
			}
			double xL = sys->nodeInfoVecs.nodeLocX[idA];
			double yL = sys->nodeInfoVecs.nodeLocY[idA];
			double zL = sys->nodeInfoVecs.nodeLocZ[idA];
			double xR = sys->nodeInfoVecs.nodeLocX[idB];
			double yR = sys->nodeInfoVecs.nodeLocY[idB];
			double zR = sys->nodeInfoVecs.nodeLocZ[idB];

			double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
			double strain = (L1 - L0) / L0;
			ofs << std::fixed << strain   << std::endl;

		}
		ofs.close();

	}

	//now print platelets 
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation+"/Platelet";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->generalParams.maxPltCount;
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " << maxPltCount  << " float" << std::endl;
		for (unsigned i = 0; i< maxPltCount; i++) {
			xPos = sys->pltInfoVecs.pltLocX[i];
			yPos = sys->pltInfoVecs.pltLocY[i];
			zPos = sys->pltInfoVecs.pltLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;

		unsigned numNumsInCells = 1 + maxPltCount;

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		
		//place edges as cells of type 2. 
		ofs<< maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points
				
		ofs << 2 << std::endl;//scatter points for capsid
		
	}
	
	//now print platelets with attatchments
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation+"/PlateletConn";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->pltInfoVecs.pltLocX.size();

		unsigned num_connections = sys->pltInfoVecs.numConnections;
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		
		ofs << "POINTS " << maxPltCount + num_connections << " float" << std::endl;
		for (unsigned i = 0; i< maxPltCount; i++) {
			xPos = sys->pltInfoVecs.pltLocX[i];
			yPos = sys->pltInfoVecs.pltLocY[i];
			zPos = sys->pltInfoVecs.pltLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//set location for nodes that plt is connected to
		//ie  
		for (unsigned i = 0; i < num_connections; i++ ) {
			unsigned node_id = sys->pltInfoVecs.nodeImagingConnection[i];
			if (node_id < sys->generalParams.maxNodeCount) {
				xPos = sys->nodeInfoVecs.nodeLocX[node_id];
				yPos = sys->nodeInfoVecs.nodeLocY[node_id];
				zPos = sys->nodeInfoVecs.nodeLocZ[node_id];
				
				ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
			}
			else{
				std::cout<<"imaging not working, node out of bounds" << std::endl;
			}
		}


		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;
		numCells += num_connections;//add conections cells for edges

		unsigned numNumsInCells = 1 + maxPltCount;
		numNumsInCells += 3 * num_connections;//3 numbers per edge

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		
		//place edges as cells of type 2. 
		ofs<< maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		
		for (unsigned edge = 0; edge < num_connections; edge++ ){

			//because nodes are placed after platelets, their vtk file id is incremented. 
			//notice that this represents the vtk id, not the id within the c++ program
			unsigned node_id_vtk = maxPltCount + edge;
			unsigned node_id = sys->pltInfoVecs.nodeImagingConnection[edge];
			
			unsigned plt_id = sys->pltInfoVecs.pltImagingConnection[edge];
				
			ofs <<2<< " "<< node_id_vtk << " "<< plt_id <<std::endl;

		}
		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points
				
		ofs << 2 << std::endl;//scatter points for capsid
		
	//	std::cout<<'here3'<<std::flush;
		for (unsigned edge = 0; edge< num_connections; edge++ ){
			ofs<< 3 <<std::endl;
		}
		ofs.close();
	}

};
