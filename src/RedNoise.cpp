#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <limits>

#define WIDTH 320
#define HEIGHT 240

const float minT = 2.60102e-12;

// #define WIDTH 512
// #define HEIGHT 512

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfSteps) {
	std::vector<float> result;
	float ratio = (to - from) / (numberOfSteps - 1);
	result.push_back(from);
	for (size_t i = 1; i < numberOfSteps; i++) {
		result.push_back(result[i - 1] + ratio);
	}
	return result;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	// 	float numberOfSteps = std::max(std::max(abs(xDiff), abs(yDiff)), abs(zDiff)); ??? --- ???
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	float zStepSize = zDiff/numberOfSteps;

	// std::vector<float> depthPoints = interpolateDepthPoints(from, to, numberOfSteps)
	for (float i=0.0; i<numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		float z = -(from.depth + (zStepSize*i));

		uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		if (z == 0) {
			// Can be deleted.
			window.setPixelColour(round(x), round(y), RGBColour);
		}
		else if (1.0/z > depthBuffer[round(x)][round(y)]) {
			window.setPixelColour(round(x), round(y), RGBColour);
			depthBuffer[round(x)][round(y)] = 1.0/z;
		}
		
	}
}

void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	drawLine(window, triangle.v0(), triangle.v1(), colour, depthBuffer);
	drawLine(window, triangle.v1(), triangle.v2(), colour, depthBuffer);
	drawLine(window, triangle.v2(), triangle.v0(), colour, depthBuffer);
}

void sortTriangleVertices(std::array<CanvasPoint, 3> &vertices) {
	if (vertices[1].y < vertices[0].y) {
		std::swap(vertices[1], vertices[0]);
	}
	if (vertices[2].y < vertices[1].y) {
		std::swap(vertices[2], vertices[1]);
	}
	if (vertices[1].y < vertices[0].y) {
		std:: swap(vertices[1], vertices[0]);
	}
}

void drawFlatBottomTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	float m0 = (v1.x - v0.x) / (v1.y - v0.y);
	float m1 = (v2.x - v0.x) / (v2.y - v0.y);

	float m2 = (v1.depth - v0.depth) / (v1.y - v0.y);
	float m3 = (v2.depth - v0.depth) / (v2.y - v0.y);

	// int yStart = round(v0.y);
	// int yEnd = round(v2.y);

	for (float y = v0.y; y < v2.y; y++) {
		float x0 = m0 * ((float)y - v0.y) + v0.x;
		float x1 = m1 * ((float)y - v0.y) + v0.x;

		float z0 = m2 * ((float)y - v0.y) + v0.depth;
		float z1 = m3 * ((float)y - v0.y) + v0.depth;

		drawLine(window, CanvasPoint(x0, y, z0), CanvasPoint(x1, y, z1), colour, depthBuffer);

		// int xStart = round(x0);
		// int xEnd = round(x1);

		// for (int x = xStart; x < xEnd; x++) {
		// 	uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		// 	window.setPixelColour(x, y, RGBColour);		
		// }
	}
}

void drawFlatTopTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	float m0 = (v2.x - v0.x) / (v2.y - v0.y);
	float m1 = (v2.x - v1.x) / (v2.y - v1.y);
	
	float m2 = (v2.depth - v0.depth) / (v2.y - v0.y);
	float m3 = (v2.depth - v1.depth) / (v2.y - v1.y);


	// int yStart = round(v0.y);
	// int yEnd = round(v2.y);

	for (float y = v0.y; y < v2.y; y++) {
		float x0 = m0 * ((float)y - v0.y) + v0.x;
		float x1 = m1 * ((float)y - v1.y) + v1.x;

		float z0 = m2 * ((float)y - v0.y) + v0.depth;
		float z1 = m3 * ((float)y - v1.y) + v1.depth;

		drawLine(window, CanvasPoint(x0, y, z0), CanvasPoint(x1, y, z1), colour, depthBuffer);


		// int xStart = round(x0);
		// int xEnd = round(x1);

		// for (int x = xStart; x < xEnd; x++) {
		// 	uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		// 	window.setPixelColour(x, y, RGBColour);
		// }
	}
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	drawStrokedTriangle(window, triangle, colour, depthBuffer);

	// Sort vertices by y-coordinate.
	sortTriangleVertices(triangle.vertices);

	if (triangle.vertices[0].y == triangle.vertices[1].y) { // flat top triangle
		if (triangle.v1().x < triangle.v0().x) {
			std::swap(triangle.v1(), triangle.v0());
		}
		drawFlatTopTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), colour, depthBuffer);
	}
	else if (triangle.vertices[1].y == triangle.vertices[2].y) { // flat bottom triangle
		if (triangle.v2().x < triangle.v1().x) {
			std::swap(triangle.v2(), triangle.v1());
		}
		drawFlatBottomTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), colour, depthBuffer);
	}
	else { // general triangle
		float ratio = (triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y);
		
		CanvasPoint vi(0.0, 0.0, 0.0);
		vi.x = triangle.v0().x + (triangle.v2().x - triangle.v0().x) * ratio;
		vi.y = triangle.v0().y + (triangle.v2().y - triangle.v0().y) * ratio;
		vi.depth = triangle.v0().depth + (triangle.v2().depth - triangle.v0().depth) * ratio;

		if (triangle.v1().x < vi.x) { // vi on the right side
			drawFlatBottomTriangle(window, triangle.v0(), triangle.v1(), vi, colour, depthBuffer);
			drawFlatTopTriangle(window, triangle.v1(), vi, triangle.v2(), colour, depthBuffer);
		}
		else { // vi on the left side
			drawFlatBottomTriangle(window, triangle.v0(), vi, triangle.v1(), colour, depthBuffer);
			drawFlatTopTriangle(window, vi, triangle.v1(), triangle.v2(), colour, depthBuffer);
		}
	}
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			// float red = rand() % 256;
			float red = 0.0;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

int getColourIndexByName(std::vector<Colour> colourPalette, std::string colourName) {
	for (int i = 0; i < colourPalette.size(); i++) {
		if (colourPalette[i].name == colourName) {
			return i;
		}
	}
	return -1;
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, bool &enableOrbit) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_o) {
			enableOrbit = !enableOrbit;
		}
		else if (event.key.keysym.sym == SDLK_LEFT) {
			cameraPosition.x -= 1;
			std::cout << "LEFT" << std::endl;
			std::cout << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			cameraPosition.x += 1;
			std::cout << "RIGHT" << std::endl;
			std::cout << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;

		}
		else if (event.key.keysym.sym == SDLK_UP) {
			cameraPosition.z -= 1;
			std::cout << "UP" << std::endl;
			std::cout << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;

		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			cameraPosition.z += 1;
			std::cout << "DOWN" << std::endl;
			std::cout << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
		}			
		else if (event.key.keysym.sym == SDLK_u) {
			CanvasPoint v0(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v1(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint v2(rand() % WIDTH, rand() % HEIGHT);
			
			CanvasTriangle triangle(v0, v1, v2);

			Colour colour(rand() % 256, rand() % 256, rand() % 256);
			std::vector<std::vector<float> > depthBuffer(WIDTH, std::vector<float> (HEIGHT, 0));

			drawFilledTriangle(window, triangle, colour, depthBuffer);

			std::cout << "DRAW FILLED TRIANGLE" << std::endl;
		}
	} 
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength) {
	CanvasPoint canvasIntersectionPoint;
	
	vertexPosition = vertexPosition - cameraPosition;

	float u =  100.0 * (-1) * focalLength * (vertexPosition.x / vertexPosition.z) + (WIDTH / 2);
	float v =  100.0 * focalLength * (vertexPosition.y / vertexPosition.z) + (HEIGHT / 2);

	canvasIntersectionPoint.x = u;
	canvasIntersectionPoint.y = v;
	canvasIntersectionPoint.depth = vertexPosition.z;

	return canvasIntersectionPoint;
}

void parseOBJFile(const std::string &filename, std::vector<Colour> colourPalette, std::vector<ModelTriangle> &modelTriangles) {
	const float scalingFactor = 0.17; 

	std::vector<glm::vec3> vertices;
	std::ifstream inputStream(filename);
	std::string nextLine;

	Colour colour;

	int indexObject = 0;

	// Get the "mtllib cornell-box.mtl"
	std::getline(inputStream, nextLine);

	// Get the empty line
	std::getline(inputStream, nextLine);

	while(std::getline(inputStream, nextLine)) {
		// std::cout << nextLine << std::endl;

		if (nextLine[0] == 'o') {

		}
		else if (nextLine[0] == 'u' && nextLine[1] == 's') { // bla bla ... change this!
			std::vector<std::string> usemtlColour = split(nextLine, ' ');
			colour = colourPalette[getColourIndexByName(colourPalette, usemtlColour[1])];
		}
		else if (nextLine[0] == 'v') {
			std::vector<std::string> vertex = split(nextLine, ' ');
			vertices.push_back(scalingFactor * glm::vec3(std::stof(vertex[1]), std::stof(vertex[2]), std::stof(vertex[3])));
		}
		else if (nextLine[0] == 'f') {
			std::vector<std::string> facets = split(nextLine, ' ');
			// facets = split(nextLine, '/');

			// std::cout << facets[0] << " " << facets[1].substr(0, facets[1].size() - 1) << " " << facets[2] << " " << facets[3] << std::endl;

			ModelTriangle triangle;

			triangle.vertices[0] = vertices[std::stoi(facets[1].substr(0, facets[1].size() - 1)) - 1];
			triangle.vertices[1] = vertices[std::stoi(facets[2].substr(0, facets[2].size() - 1)) - 1];
			triangle.vertices[2] = vertices[std::stoi(facets[3].substr(0, facets[3].size() - 1)) - 1];

			triangle.colour = colour;

			modelTriangles.push_back(triangle);
		}
		else {
			indexObject++; // TODO: parse in a similiar way to .mtl files
		}
	}
}

void readMTLFile(const std::string &filename, std::vector<Colour> &colourPalette) {
	std::ifstream inputStream(filename);
	std::string nextLine;
	Colour colour;

	while(std::getline(inputStream, nextLine)) {

		std::vector<std::string> line = split(nextLine, ' ');

		if (line[0] == "newmtl") {
			colour.name = line[1];
		}
		else if (line[0] == "Kd") {
			colour.red = 255 * std::stof(line[1]);
			colour.green = 255 * std::stof(line[2]);
			colour.blue = 255 * std::stof(line[3]);

			if (inputStream.eof()) {
				colourPalette.push_back(colour);
				return;
			}
		}
		else {
			colourPalette.push_back(colour);
		}
	}
}

void clearDepthBuffer(std::vector<std::vector<float> > &depthBuffer) {
	for (int i = 0; i < depthBuffer.size(); i++) {
		std::fill(depthBuffer[i].begin(), depthBuffer[i].end(), 0);
	}
}

glm::mat3 lookAt(glm::vec3 from, glm::vec3 to) {
	glm::vec3 forward = normalize(from - to);
	glm::vec3 right = cross(glm::vec3(0.0, 1.0, 0.0), forward);
	glm::vec3 up = cross(forward, right);

	glm::mat3 cameraOrientation = glm::mat3(right, up, forward);

	return cameraOrientation;
}

void drawRasterizedScene(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, std::vector<std::vector<float> > &depthBuffer) {
	window.clearPixels();
	clearDepthBuffer(depthBuffer);

	CanvasPoint v0;
	CanvasPoint v1;
	CanvasPoint v2;

	for (int i = 0; i < modelTriangles.size(); i++) {
		for (int j = 0; j < modelTriangles[i].vertices.size(); j++) {
			glm::vec3 cameraToVertex = modelTriangles[i].vertices[j] - cameraPosition;
			glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;
			modelTriangles[i].vertices[j] = adjustedVector + cameraPosition; // - cameraPosition
		}

		v0 = getCanvasIntersectionPoint(cameraPosition, modelTriangles[i].vertices[0], 2);
		v1 = getCanvasIntersectionPoint(cameraPosition, modelTriangles[i].vertices[1], 2);
		v2 = getCanvasIntersectionPoint(cameraPosition, modelTriangles[i].vertices[2], 2);

		CanvasTriangle triangle(v0, v1, v2);
		drawFilledTriangle(window, triangle, modelTriangles[i].colour, depthBuffer);
	}
}

bool validIntersection(float t, float u, float v) {
	return (t > 0.0) && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0);
}

glm::vec3 getIntersectionPoint(ModelTriangle triangle, float u, float v) {
	glm::vec3 intersectionPoint;
	intersectionPoint = triangle.vertices[0] + u * (triangle.vertices[1] - triangle.vertices[0]) + v * (triangle.vertices[2] - triangle.vertices[0]);
	return intersectionPoint; // Maybe origin + intersectionPoint???
}

RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> modelTriangles, glm::mat3 cameraOrientation, bool shadowRay) {
	RayTriangleIntersection rayTriangleIntersection;

	rayTriangleIntersection.distanceFromCamera = std::numeric_limits<float>::max();


	for (int i = 0; i < modelTriangles.size(); i++) {
		
		// adjust triangle to rotation perspective - for orbiting
		ModelTriangle triangle = modelTriangles[i];

		// if (! shadowRay) {
		// 	for (int j = 0; j < triangle.vertices.size(); j++) {
		// 		glm::vec3 cameraToVertex = triangle.vertices[j] - cameraPosition;
		// 		glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;
		// 		triangle.vertices[j] = adjustedVector + cameraPosition; // - cameraPosition
		// 	}
		// }

		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;

		if (validIntersection(t, u, v) && t > 0.002) { // Change t > 0.01 ! // 0.0018
			if (t < rayTriangleIntersection.distanceFromCamera) {
				rayTriangleIntersection.distanceFromCamera = t;
				rayTriangleIntersection.intersectedTriangle = triangle;
				rayTriangleIntersection.intersectionPoint = getIntersectionPoint(triangle, u, v);
				rayTriangleIntersection.triangleIndex = i;
			}
		}
	}
	return rayTriangleIntersection;
}

void drawRayTracedScene(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, glm::vec3 cameraPosition, glm::vec3 rayDirection, glm::mat3 cameraOrientation, glm::vec3 lightPosition) {
	window.clearPixels();

	float imageAspectRatio = window.width / (float)window.height; // assuming width > height 
	// float minDistance = std::numeric_limits<float>::max();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			rayDirection.x = (2 * ((x + 0.5) / window.width) - 1) * imageAspectRatio; 
			rayDirection.y = (1 - 2 * ((y + 0.5) / window.height)); 

			// std::cout << rayDirection.x << " " << rayDirection.y << " " << rayDirection.z << std::endl;

			// Make it look bigger
			rayDirection.x /= 3.0;
			rayDirection.y /= 3.0;

			rayDirection = glm::normalize(rayDirection); // Should I? TODO: Problem when zooming in and out!
			RayTriangleIntersection rayTriangleIntersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles, cameraOrientation, false);
			if (rayTriangleIntersection.distanceFromCamera != std::numeric_limits<float>::max()) {
				glm::vec3 newCameraPosition = rayTriangleIntersection.intersectionPoint;
				glm::mat3 newCameraOrientation = lookAt(newCameraPosition, lightPosition);
				glm::vec3 newRayDirection = lightPosition - newCameraPosition; // ?
				glm::normalize(newRayDirection); // Should I? TODO: Problem when zooming in and out!

				float distanceToLightPosition = glm::distance(lightPosition, newCameraPosition);
				// if (distanceToLightPosition < minDistance) {
				// 	std::cout << "d -> l = " << distanceToLightPosition << std::endl;
				// 	minDistance = distanceToLightPosition;
				// }
				RayTriangleIntersection newRayTriangleIntersection = getClosestIntersection(newCameraPosition, newRayDirection, modelTriangles, newCameraOrientation, true);

				// Light can't be seen
				if (newRayTriangleIntersection.distanceFromCamera < distanceToLightPosition && rayTriangleIntersection.triangleIndex != newRayTriangleIntersection.triangleIndex) { 
					// window.setPixelColour(x, y, 0);
				} 
				else {
					Colour RGBColour = rayTriangleIntersection.intersectedTriangle.colour;
					uint32_t colour = (255 << 24) + (int(RGBColour.red) << 16) + (int(RGBColour.green) << 8) + int(RGBColour.blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
	// std::cout << "minDistance=" << minDistance << std::endl;
}

// Rotate camera about X axis with a given rotation angle
glm::vec3 rotateCameraAboutXAxis(glm::vec3 cameraPosition, float angle) {
	glm::mat3 rotationMatrix = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, cos(angle), sin(angle),
		0.0, -sin(angle), cos(angle)
	);
	glm::vec3 newCameraPosition = rotationMatrix * cameraPosition;
	return newCameraPosition;
}

// Rotate camera about Y axis with a given rotation angle
glm::vec3 rotateCameraAboutYAxis(glm::vec3 cameraPosition, float angle) {
	glm::mat3 rotationMatrix = glm::mat3(
		cos(angle), 0.0, -sin(angle),
		0.0, 1.0, 0,
		sin(angle), 0.0, cos(angle)
	);
	glm::vec3 newCameraPosition = rotationMatrix * cameraPosition;
	return newCameraPosition;
}


int main(int argc, char *argv[]) {
	std::string objFilename = "./models/cornell-box.obj";
	std::string mtlFilename = "./models/cornell-box.mtl";

	std::vector<Colour> colourPalette;
	std::vector<ModelTriangle> modelTriangles;
	std::vector<std::vector<float> > depthBuffer(WIDTH, std::vector<float> (HEIGHT, 0));

	readMTLFile(mtlFilename, colourPalette);
	parseOBJFile(objFilename, colourPalette, modelTriangles);

	glm::vec3 cameraPosition(0.0, 0.0, 2.0); // Try (0.0, 0.0, 4.0) for rasterising
	glm::vec3 rayDirection(0.0, 0.0, -1.0);
	glm::vec3 lightPosition(0.0, 0.35, 0.0);

	glm::vec3 right(1.0, 0.0, 0.0); 
	glm::vec3 up(0.0, 1.0, 0.0);
	glm::vec3 forward(0.0, 0.0, 1.0);

	glm::mat3 cameraOrientation(right, up, forward);

	// glm::vec3 column0 = glm::vec3(1.0, 0.0, 0.0);
	// glm::vec3 column1 = glm::vec3(0.0, 0.7, 0.7);
	// glm::vec3 column2 = glm::vec3(0.0, -0.7, 0.7);
	// glm::mat3 n = glm::mat3(column0, column1, column2); 

	glm::mat3 rotationDefault = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0, 0.0
	);

	// cameraOrientation = rotationDefault * cameraOrientation;

	// std::cout << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z;

	// for (int i = 0; i < modelTriangles.size(); i++) {
	// 	for (int j = 0; j < modelTriangles[i].vertices.size(); j++) {
	// 		glm::vec3 cameraToVertex = modelTriangles[i].vertices[j] - cameraPosition;
	// 		glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;
	// 		modelTriangles[i].vertices[j] = adjustedVector + cameraPosition;
	// 	}
	// }

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	// CanvasPoint from(0, 0);
	// CanvasPoint to(WIDTH, HEIGHT);

	// CanvasPoint v0(30, 30);
	// CanvasPoint v1(30, 60);
	// CanvasPoint v2(60, 60);
	// CanvasTriangle triangle(v0, v1, v2);

	// CanvasPoint v0Test(120, 180);
	// CanvasPoint v1Test(20, 10);
	// CanvasPoint v2Test(5, 90);
	// CanvasTriangle triangleTest(v0Test, v1Test, v2Test);

	// CanvasPoint v00(164.415, 101.147, -3.02656);
	// CanvasPoint v11(98.4548, 102.272, -3.21356);
	// CanvasPoint v22(109.784, 100.6, -2.94326);
	// CanvasTriangle triangleTest00(v00, v11, v22);

	// Colour BLACK(255, 255, 255);

	bool enableOrbit = false;
	bool upDown = true;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cameraPosition, cameraOrientation, enableOrbit);
		drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
		// drawRasterizedScene(window, modelTriangles, cameraPosition, cameraOrientation, depthBuffer);

		// if (upDown) {
		// 	lightPosition.y += 0.005;
		// 	if (lightPosition.y > 0.46) {
		// 		upDown = false;
		// 	}
		// }
		// else {
		// 	lightPosition.y -= 0.005;
		// 	if (lightPosition.y < 0.0) {
		// 		upDown = true;
		// 	}
		// }

		if (enableOrbit) {
			// float angle = 3.14 / 2048.0;
			float angle = 3.14 / 16.0;
			cameraPosition = rotateCameraAboutYAxis(cameraPosition, angle);
			// cameraPosition = rotateCameraAboutXAxis(cameraPosition, angle);
			cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
		}
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
