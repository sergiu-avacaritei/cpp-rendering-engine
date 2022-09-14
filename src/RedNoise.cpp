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
#include <TextureMap.h>
#include <TexturePoint.h>

// #define WIDTH 320
// #define HEIGHT 240

#define WIDTH 480
#define HEIGHT 360

enum RenderOption { Wireframe, Rasterise, Raytrace } renderOption;

bool orbiting = false;

TextureMap textureMap("./models/texture-logo.ppm");
TextureMap chessTextureMap("./models/check3.ppm");
TextureMap cobblesTextureMap("./models/texture.ppm");


// const float minT = 2.60102e-12;

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, std::vector<std::vector<float> > &depthBuffer) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	
	float numberOfSteps = std::max(fabs(xDiff), fabs(yDiff));
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

int getColourIndexByName(std::vector<Colour> colourPalette, std::string colourName) {
	for (int i = 0; i < colourPalette.size(); i++) {
		if (colourPalette[i].name == colourName) {
			return i;
		}
	}
	return -1;
}


// -------------------------------------------- TEXTURE MAPPING RASTERIZER --------------------------- //

void drawTexturedLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, std::vector<std::vector<uint32_t> > texture) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;

	float xtDiff = to.texturePoint.x - from.texturePoint.x;
	float ytDiff = to.texturePoint.y - from.texturePoint.y;

	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));

	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;

	float xtStepSize = xtDiff/numberOfSteps;
	float ytStepSize = ytDiff/numberOfSteps;

	for (float i=0.0; i<numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);

		float xt = from.texturePoint.x + (xtStepSize * i);
		float yt = from.texturePoint.y + (ytStepSize * i);
		// uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		// if (x > 0 && x < WIDTH && y > 0 && y < HEIGHT) {
			window.setPixelColour(round(x), round(y), texture[round(yt)][round(xt)]);


		// }
		// window.setPixelColour(round(x), round(y), texture[round(yt)][round(xt)]);
	}

}

void drawFlatBottomTexturedTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, std::vector<std::vector<uint32_t> > texture) {
	float m0 = (v1.x - v0.x) / (v1.y - v0.y);
	float m1 = (v2.x - v0.x) / (v2.y - v0.y);

	float m0_t = (v1.texturePoint.x - v0.texturePoint.x) / (v1.texturePoint.y - v0.texturePoint.y);
	float m1_t = (v2.texturePoint.x - v0.texturePoint.x) / (v2.texturePoint.y - v0.texturePoint.y);

	

	// std::vector<uint32_t> textureLine;

	// int yStart = round(v0.y);
	// int yEnd = round(v2.y);

	// for (int y = v0.texturePoint.y; y < v2.texturePoint.y; y++) {
	// 	float x0_t = m0 * ((float)y - v0.texturePoint.y) + v0.texturePoint.x;
	// 	float x1_t = m1 * ((float)y - v0.texturePoint.y) + v0.texturePoint.x;

	// }

	for (float y = v0.y, yt = v0.texturePoint.y; y < v2.y; y+= 1.0) {
		float x0 = m0 * ((float)y - v0.y) + v0.x;
		float x1 = m1 * ((float)y - v0.y) + v0.x;
			
		float x0_t = m0_t * ((float)yt - v0.texturePoint.y) + v0.texturePoint.x;
		float x1_t = m1_t * ((float)yt - v0.texturePoint.y) + v0.texturePoint.x;

		CanvasPoint p1, p2;

		p1.x = x0;
		p1.y = y;
		p1.texturePoint.x = x0_t;
		p1.texturePoint.y = yt;

		p2.x = x1;
		p2.y = y;
		p2.texturePoint.x = x1_t;
		p2.texturePoint.y = yt;

		// textureLine = interpolateTexturePoints(p1.texturePoint, p2.texturePoint, stepCountC);
		drawTexturedLine(window, p1, p2, texture);

		yt = yt + (v2.texturePoint.y - v0.texturePoint.y) / (v2.y - v0.y);
		// drawLine(window, p1, p2, Colour(255,0,0));

		// int xStart = round(x0);
		// int xEnd = round(x1);

		// for (int x = xStart; x < xEnd; x++) {
		// 	uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		// 	window.setPixelColour(x, y, RGBColour);		
		// }
	}
}

void drawFlatTopTexturedTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, std::vector<std::vector<uint32_t> > texture) {
	float m0 = (v2.x - v0.x) / (v2.y - v0.y);
	float m1 = (v2.x - v1.x) / (v2.y - v1.y);

	float m0_t = (v2.texturePoint.x - v0.texturePoint.x) / (v2.texturePoint.y - v0.texturePoint.y);
	float m1_t = (v2.texturePoint.x - v1.texturePoint.x) / (v2.texturePoint.y - v1.texturePoint.y);

	// int yStart = round(v0.y);
	// int yEnd = round(v2.y);

	// std::cout << "ratio = " << (v2.texturePoint.y - v0.texturePoint.y) / (v2.y - v0.y);

	for (float y = v0.y, yt = v0.texturePoint.y; y < v2.y; y+= 1.0) {
		float x0 = m0 * ((float)y - v0.y) + v0.x;
		float x1 = m1 * ((float)y - v1.y) + v1.x;

		float x0_t = m0_t * ((float)yt - v0.texturePoint.y) + v0.texturePoint.x;
		float x1_t = m1_t * ((float)yt - v1.texturePoint.y) + v1.texturePoint.x;

		CanvasPoint p1, p2;

		p1.x = x0;
		p1.y = y;
		p1.texturePoint.x = x0_t;
		p1.texturePoint.y = yt;

		p2.x = x1;
		p2.y = y;
		p2.texturePoint.x = x1_t;
		p2.texturePoint.y = yt;

		drawTexturedLine(window, p1, p2, texture);

		yt = yt + (v2.texturePoint.y - v0.texturePoint.y) / (v2.y - v0.y);
		// drawLine(window, p1, p2, Colour(255,0,0));
		// uint32_t RGBColour = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(255);
		// window.setPixelColour(v1.x, v1.y, RGBColour);



		// int xStart = round(x0);
		// int xEnd = round(x1);

		// for (int x = xStart; x < xEnd; x++) {
		// 	uint32_t RGBColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		// 	window.setPixelColour(x, y, RGBColour);
		// }
	}

}

void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, CanvasTriangle texturedTriangle) {
	// Load texture.ppm
	// TextureMap textureMap("texture.ppm");

	// 1D -> 2D texture
	std::vector<std::vector<uint32_t> > texture(cobblesTextureMap.width, std::vector<uint32_t> (cobblesTextureMap.height));
	for (int y = 0; y < cobblesTextureMap.height; y++) {
		for (int x = 0; x < cobblesTextureMap.width; x++) {
			texture[y][x] = cobblesTextureMap.pixels[y * cobblesTextureMap.width + x];
		}
	}

	// sortTriangleVertices(triangle.vertices);
	// sortTriangleVertices(texturedTriangle.vertices);


	// Populate texturePoint field of every vertex
	triangle.v0().texturePoint.x = texturedTriangle.v0().x;
	triangle.v1().texturePoint.x = texturedTriangle.v1().x;
	triangle.v2().texturePoint.x = texturedTriangle.v2().x;

	triangle.v0().texturePoint.y = texturedTriangle.v0().y;
	triangle.v1().texturePoint.y = texturedTriangle.v1().y;
	triangle.v2().texturePoint.y = texturedTriangle.v2().y;

	// sortTriangleVertices(triangle.vertices);

	// Sort vertices
	sortTriangleVertices(triangle.vertices);
	// sortTriangleVertices(texturedTriangle.vertices);

	float ratio = (triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y);
	// float t_ratio = (texturedTriangle.v1().y - texturedTriangle.v0().y) / (texturedTriangle.v2().y - texturedTriangle.v0().y);
		
	CanvasPoint vi(0, 0);
	vi.x = triangle.v0().x + (triangle.v2().x - triangle.v0().x) * ratio;
	vi.y = triangle.v0().y + (triangle.v2().y - triangle.v0().y) * ratio;

	// sortTriangleVertices(texturedTriangle.vertices); // Should I sort them? !!!
	// vi.texturePoint.x = texturedTriangle.v0().x + (texturedTriangle.v2().x - texturedTriangle.v0().x) * ratio;
	// vi.texturePoint.y = texturedTriangle.v0().y + (texturedTriangle.v2().y - texturedTriangle.v0().y) * ratio;

	// std::cout << triangle.v0().texturePoint << " " << vi.texturePoint << std::endl;

	vi.texturePoint.x = triangle.v0().texturePoint.x + (triangle.v2().texturePoint.x - triangle.v0().texturePoint.x) * ratio;
	vi.texturePoint.y = triangle.v0().texturePoint.y + (triangle.v2().texturePoint.y - triangle.v0().texturePoint.y) * ratio;

	// std::cout << vi.texturePoint.x << " " << vi.texturePoint.y << std::endl;
	// drawStrokedTriangle(window, triangle, Colour(255,255,255));

	// std::cout << "leftPoint.y = " << triangle.v1().texturePoint.y << " || " << "intersectionPoint.y = " << vi.texturePoint.y << " || " << "bottomPoint.y = " << triangle.v2().texturePoint.y << std::endl;

	if (triangle.v1().x < vi.x) { // vi on the right side
			drawFlatBottomTexturedTriangle(window, triangle.v0(), triangle.v1(), vi, texture);
			drawFlatTopTexturedTriangle(window, triangle.v1(), vi, triangle.v2(), texture);
	}
	else { // vi on the left side
			drawFlatBottomTexturedTriangle(window, triangle.v0(), vi, triangle.v1(), texture);
			drawFlatTopTexturedTriangle(window, vi, triangle.v1(), triangle.v2(), texture);
	}

}

// -------------------------------------------- TEXTURE MAPPING RASTERIZER --------------------------- //

glm::mat3 lookAt(glm::vec3 from, glm::vec3 to) {
	glm::vec3 forward = normalize(from - to);
	glm::vec3 right = cross(glm::vec3(0.0, 1.0, 0.0), forward);
	glm::vec3 up = cross(forward, right);

	glm::mat3 cameraOrientation = glm::mat3(right, up, forward);

	return cameraOrientation;
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

glm::vec3 calculateTriangleNormal(ModelTriangle triangle) {
	glm::vec3 triangleNormal = glm::cross(triangle.vertices[1] - triangle.vertices[0], triangle.vertices[2] - triangle.vertices[0]);
	return triangleNormal;
}

std::string objectType = "box";

void parseOBJFile(const std::string &filename, std::vector<Colour> colourPalette, std::vector<ModelTriangle> &modelTriangles) {
	const float scalingFactor = 0.17; 
		// const float scalingFactor = 0.7; for loggo


	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> verticesNormal;
	std::vector<TexturePoint> textureVertices;
	std::ifstream inputStream(filename);
	std::string nextLine;

	Colour colour(255, 0, 0);

	int indexObject = 0;

	// Get the "mtllib cornell-box.mtl"
	// std::getline(inputStream, nextLine);

	// Get the empty line
	// std::getline(inputStream, nextLine);

	while(std::getline(inputStream, nextLine)) {
		// std::cout << nextLine << std::endl;

		if (nextLine[0] == 'o') {
			std::vector<std::string> object = split(nextLine, ' ');
			
			if (object[1] == "sphere") {
				objectType = object[1];
			}
			std::cout << objectType << std::endl;
		}
		else if (nextLine[0] == 'u' && nextLine[1] == 's') { // bla bla ... change this!
			std::vector<std::string> usemtlColour = split(nextLine, ' ');
			colour = colourPalette[getColourIndexByName(colourPalette, usemtlColour[1])];

		}
		else if (nextLine[0] == 'v' && nextLine[1] != 'n' && nextLine[1] != 't') {
			// glm::vec3 origin(0.0, 0.0, 0.0);
			std::vector<std::string> vertex = split(nextLine, ' ');
			
			// TODO: CHANGE IT FOR CORNELL BOX
			// glm::vec3 origin(0.0,-1.5,0.0); // This is for sphere.obj
			glm::vec3 origin(0.0,0.0,0.0);
			vertices.push_back(scalingFactor * (origin + glm::vec3(std::stof(vertex[1]), std::stof(vertex[2]), std::stof(vertex[3])))); // /= 505.0f for logo
		}
		else if (nextLine[0] == 'v' && nextLine[1] == 't') { // vt
			std::vector<std::string> vertexTexture = split(nextLine, ' ');
			textureVertices.push_back(TexturePoint(std::stof(vertexTexture[1]), std::stof(vertexTexture[2])));
		}
		else if (nextLine[0] == 'f') {
			std::vector<std::string> facets = split(nextLine, ' ');

			std::vector<std::string> f0 = split(facets[1], '/');
			std::vector<std::string> f1 = split(facets[2], '/');
			std::vector<std::string> f2 = split(facets[3], '/');

			// std::cout << f0[0] << std::endl;

			// facets = split(nextLine, '/');
			
			// std::cout << facets[1].size() << std::endl;
			if (facets[1][facets[1].size() / 2] == '/' && facets[1][facets[1].size() / 2 - 1] == '/') {
				facets[1] = facets[1].substr(0, facets[1].size() / 2);
				facets[2] = facets[2].substr(0, facets[2].size() / 2);
				facets[3] = facets[3].substr(0, facets[3].size() / 2);
				// std::cout << "OK" << std::endl;
				// std::cout << facets[1].size() << std::endl;
			}

			// std::cout << facets[0] << " " << facets[1]<< " " << facets[2] << " " << facets[3] << std::endl;

			ModelTriangle triangle;

			triangle.vertices[0] = vertices[std::stoi(f0[0]) - 1];
			triangle.vertices[1] = vertices[std::stoi(f1[0]) - 1];
			triangle.vertices[2] = vertices[std::stoi(f2[0]) - 1];

			if (colour.name == "Logo") {
				triangle.texturePoints[0] = TexturePoint(textureMap.width * textureVertices[std::stoi(f0[1]) - 1].x, textureMap.height - textureMap.height * textureVertices[std::stoi(f0[1]) - 1].y);
				triangle.texturePoints[1] = TexturePoint(textureMap.width * textureVertices[std::stoi(f1[1]) - 1].x, textureMap.height - textureMap.height * textureVertices[std::stoi(f1[1]) - 1].y);
				triangle.texturePoints[2] = TexturePoint(textureMap.width * textureVertices[std::stoi(f2[1]) - 1].x, textureMap.height - textureMap.height * textureVertices[std::stoi(f2[1]) - 1].y);
				triangle.hasTexture = true;
			}

			if (colour.name == "Chess") {
				triangle.texturePoints[0] = TexturePoint(chessTextureMap.width * textureVertices[std::stoi(f0[1]) - 1].x, chessTextureMap.height - chessTextureMap.height * textureVertices[std::stoi(f0[1]) - 1].y);
				triangle.texturePoints[1] = TexturePoint(chessTextureMap.width * textureVertices[std::stoi(f1[1]) - 1].x, chessTextureMap.height - chessTextureMap.height * textureVertices[std::stoi(f1[1]) - 1].y);
				triangle.texturePoints[2] = TexturePoint(chessTextureMap.width * textureVertices[std::stoi(f2[1]) - 1].x, chessTextureMap.height - chessTextureMap.height * textureVertices[std::stoi(f2[1]) - 1].y);
				triangle.hasTexture = true;
			}

			if (colour.name == "Cobbles") {
				triangle.texturePoints[0] = TexturePoint(cobblesTextureMap.width * textureVertices[std::stoi(f0[1]) - 1].x, cobblesTextureMap.height - cobblesTextureMap.height * textureVertices[std::stoi(f0[1]) - 1].y);
				triangle.texturePoints[1] = TexturePoint(cobblesTextureMap.width * textureVertices[std::stoi(f1[1]) - 1].x, cobblesTextureMap.height - cobblesTextureMap.height * textureVertices[std::stoi(f1[1]) - 1].y);
				triangle.texturePoints[2] = TexturePoint(cobblesTextureMap.width * textureVertices[std::stoi(f2[1]) - 1].x, cobblesTextureMap.height - cobblesTextureMap.height * textureVertices[std::stoi(f2[1]) - 1].y);
				triangle.hasTexture = true;
			}

			if (colour.name == "Mirror") {
				triangle.isReflective = true;
			}

			if (objectType == "sphere") {
				triangle.belongsToSphere = true;
			}

			// std::cout << facets[1].substr(0, facets[1].size() - 1) << " " << facets[2].substr(0, facets[2].size() - 1) << " " << facets[3].substr(0, facets[3].size() - 1) << std::endl;

			triangle.normal = calculateTriangleNormal(triangle);
			triangle.normal = normalize(triangle.normal);
			
			// if (verticesNormal.size() > 0) {
			// 	std::cout << "HERE!";
			// 	triangle.hasVertexNormals = true;

			// 	triangle.verticesNormal[0] = verticesNormal[std::stoi(facets[1].substr(0, facets[1].size() - 1)) - 1];
			// 	triangle.verticesNormal[1] = verticesNormal[std::stoi(facets[2].substr(0, facets[2].size() - 1)) - 1];
			// 	triangle.verticesNormal[2] = verticesNormal[std::stoi(facets[3].substr(0, facets[3].size() - 1)) - 1];
			// }

			triangle.colour = colour;

			if (colour.name != "Black") {
				modelTriangles.push_back(triangle);
			}
			
		}
		else {
			indexObject++; // TODO: parse in a similiar way to .mtl files
		}
	}

	// Get the vertex normals for each triangle
	for (int i = 0; i < modelTriangles.size(); i++) {
		// For each vertex vi
		for (int vi = 0; vi < modelTriangles[i].vertices.size(); vi++) {
			int vertexCount = 0;
			glm::vec3 vertex_i = modelTriangles[i].vertices[vi];
			glm::vec3 sumNormal(0.0, 0.0, 0.0);
			for (int j = 0; j < modelTriangles.size(); j++) {
				for (int vj = 0; vj < modelTriangles[j].vertices.size(); vj++) {
					glm::vec3 vertex_j = modelTriangles[j].vertices[vj];

					bool shareSameVertex = vertex_i == vertex_j;

					if(shareSameVertex) {
						sumNormal = sumNormal + modelTriangles[j].normal;
						vertexCount++;
					}
				}
			}
			modelTriangles[i].verticesNormal[vi] = normalize(sumNormal / (float) vertexCount);
		}
		modelTriangles[i].hasVertexNormals = true;
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


void drawWireframedScene(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, std::vector<std::vector<float> > &depthBuffer) {
	window.clearPixels();
	clearDepthBuffer(depthBuffer);

	CanvasPoint v0;
	CanvasPoint v1;
	CanvasPoint v2;

	// TODO make a copy of modelTriangles[i]
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
		drawStrokedTriangle(window, triangle, modelTriangles[i].colour, depthBuffer);
	}
}

void drawRasterizedScene(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, std::vector<std::vector<float> > &depthBuffer) {
	window.clearPixels();
	clearDepthBuffer(depthBuffer);

	CanvasPoint v0;
	CanvasPoint v1;
	CanvasPoint v2;

	// TODO make a copy of modelTriangles[i]
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
		if (modelTriangles[i].hasTexture) {
			CanvasTriangle texturedTriangle;
			texturedTriangle.vertices[0].x = modelTriangles[i].texturePoints[0].x;
			texturedTriangle.vertices[0].y= modelTriangles[i].texturePoints[0].y;

			texturedTriangle.vertices[1].x = modelTriangles[i].texturePoints[1].x;
			texturedTriangle.vertices[1].y= modelTriangles[i].texturePoints[1].y;

			texturedTriangle.vertices[2].x = modelTriangles[i].texturePoints[2].x;
			texturedTriangle.vertices[2].y= modelTriangles[i].texturePoints[2].y;

			drawTexturedTriangle(window, triangle, texturedTriangle);
		}
		else {
			drawFilledTriangle(window, triangle, modelTriangles[i].colour, depthBuffer);
		}
		
	}
}

bool validIntersection(float t, float u, float v) {
	return (t > 0.0) && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0);
}

glm::vec3 getIntersectionPoint(ModelTriangle triangle, float u, float v) {
	glm::vec3 intersectionPoint = triangle.vertices[0] + u * (triangle.vertices[1] - triangle.vertices[0]) + v * (triangle.vertices[2] - triangle.vertices[0]);
	return intersectionPoint; 
}

ModelTriangle adjustModelTriangleWRTCameraOrientation(ModelTriangle triangle, glm::vec3 cameraPosition, glm::mat3 cameraOrientation) {
	ModelTriangle adjustedTriangle = triangle;
	for (int i = 0; i < adjustedTriangle.vertices.size(); i++) {
		glm::vec3 cameraToVertex = adjustedTriangle.vertices[i] - cameraPosition;
		glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;
		adjustedTriangle.vertices[i] = adjustedVector + cameraPosition; 
	}
	return adjustedTriangle;
}

glm::vec3 getRayTriangleIntersection(ModelTriangle triangle, glm::vec3 cameraPosition, glm::vec3 rayDirection) {
	glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
	glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
	glm::mat3 DEMatrix(-rayDirection, e0, e1);
	glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
	return possibleSolution;
}

glm::vec3 calculateTriangleNormalFromBarycentricCoordinates(ModelTriangle triangle, glm::vec3 barycentricCoordinates) {
    glm::vec3 vertexNormal0 = glm::normalize(triangle.verticesNormal[0]);
    glm::vec3 vertexNormal1 = glm::normalize(triangle.verticesNormal[1]);
    glm::vec3 vertexNormal2 = glm::normalize(triangle.verticesNormal[2]);
    glm::vec3 triangleNormal = barycentricCoordinates.x * vertexNormal1 + barycentricCoordinates.y * vertexNormal2 + barycentricCoordinates.z * vertexNormal0;
    return triangleNormal;
}


RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> modelTriangles, glm::mat3 cameraOrientation, bool shadowRay, glm::vec3 &barycentricCoordinates) {
	RayTriangleIntersection rayTriangleIntersection;
	rayTriangleIntersection.distanceFromCamera = std::numeric_limits<float>::max();

	for (int i = 0; i < modelTriangles.size(); i++) {
		ModelTriangle triangle = modelTriangles[i];

		glm::vec3 possibleSolution = getRayTriangleIntersection(triangle, cameraPosition, rayDirection);
		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (validIntersection(t, u, v) && t > 0.002) { // Change t > 0.01 ! // 0.0018
			if (t < rayTriangleIntersection.distanceFromCamera) {
				rayTriangleIntersection.distanceFromCamera = t;
				rayTriangleIntersection.intersectedTriangle = triangle;
				rayTriangleIntersection.intersectionPoint = getIntersectionPoint(triangle, u, v);
				rayTriangleIntersection.triangleIndex = i;

				if (! shadowRay) {
					barycentricCoordinates.x = u;
					barycentricCoordinates.y = v;
					barycentricCoordinates.z = 1.0f - u - v;
				}
			}
		}
	}
	// If intersected triangle is reflective (i.e.: a mirror) then fire another ray and see what it hits
	if (rayTriangleIntersection.intersectedTriangle.isReflective && !shadowRay) {
		glm::vec3 normal = calculateTriangleNormal(rayTriangleIntersection.intersectedTriangle);
		normal = glm::normalize(normal);

		if (rayTriangleIntersection.intersectedTriangle.belongsToSphere) {
			// std::cout << "AICI" << std::endl;
			normal = calculateTriangleNormalFromBarycentricCoordinates(rayTriangleIntersection.intersectedTriangle, barycentricCoordinates);
			normal = glm::normalize(normal);
			
		}

		glm::vec3 reflectedRay = rayDirection - (2.0f * normal) * (glm::dot(rayDirection, normal));
		reflectedRay = glm::normalize(reflectedRay);
		return getClosestIntersection(rayTriangleIntersection.intersectionPoint, reflectedRay, modelTriangles, cameraOrientation, false, barycentricCoordinates);
	}
	return rayTriangleIntersection;
}

std::vector<glm::vec3> generateLightPoints(glm::vec3 lightPosition) {
	std::vector<glm::vec3> lightPoints;

	float d = 0.07f;
	float r = 0.1f;

	lightPoints.push_back(lightPosition);

	// // std::cout << "BEGIN" << std::endl;
	int count = 0;

	for (float x0 = lightPosition.x; x0 <= lightPosition.x + r; x0 = x0 + d) {
		for (float y0 = lightPosition.y; y0 <= lightPosition.y + r; y0 = y0 + d) {
			for (float z0 = lightPosition.z; z0 <= lightPosition.z + r; z0 = z0 + d) {
				if (pow(lightPosition.x - x0, 2) + pow(lightPosition.y - y0, 2) + pow(lightPosition.z - z0, 2) <= pow(r, 2)) {
					// std::cout << x0 << " " << y0 << " " << z0 << std::endl;
					lightPoints.push_back(glm::vec3(x0, y0, z0));
					count++;
				}
			}
		}
	}

	d = 0.02;

	lightPoints.push_back(glm::vec3(lightPosition.x + d, lightPosition.y, lightPosition.z));
	lightPoints.push_back(glm::vec3(lightPosition.x - d, lightPosition.y, lightPosition.z));
	lightPoints.push_back(glm::vec3(lightPosition.x, lightPosition.y + d, lightPosition.z));
	lightPoints.push_back(glm::vec3(lightPosition.x, lightPosition.y - d, lightPosition.z));
	lightPoints.push_back(glm::vec3(lightPosition.x, lightPosition.y, lightPosition.z + d));
	lightPoints.push_back(glm::vec3(lightPosition.x, lightPosition.y, lightPosition.z - d));

		
	lightPoints.push_back(glm::vec3(lightPosition.x + d, lightPosition.y + d, lightPosition.z + d));
	lightPoints.push_back(glm::vec3(lightPosition.x + d, lightPosition.y + d, lightPosition.z - d));
	lightPoints.push_back(glm::vec3(lightPosition.x + d, lightPosition.y - d, lightPosition.z + d));
	lightPoints.push_back(glm::vec3(lightPosition.x + d, lightPosition.y - d, lightPosition.z - d));
	lightPoints.push_back(glm::vec3(lightPosition.x - d, lightPosition.y - d, lightPosition.z - d));
	lightPoints.push_back(glm::vec3(lightPosition.x - d, lightPosition.y - d, lightPosition.z + d));
	lightPoints.push_back(glm::vec3(lightPosition.x - d, lightPosition.y + d, lightPosition.z - d));
	lightPoints.push_back(glm::vec3(lightPosition.x - d, lightPosition.y + d, lightPosition.z + d));

	// for (int i = 0; i < lightPoints.size(); i++) {
	// 	std::cout << lightPoints[i].x << " " << lightPoints[i].y << " " << lightPoints[i].z << std::endl;
	// }

	return lightPoints;
}

float maxDistanceToLightPosition = -1.0;
float minDistanceToLightPosition = std::numeric_limits<float>::max();

float lightRange = 3.0f; // 2.5 for cornell box

bool proximityOn = true;
bool incidenceOn = true;
bool specularOn = true;

bool gouraud = false;
bool phong = true;
bool defaultLight = false;
bool ambientLight = true;

float calculateBrightness(RayTriangleIntersection rayTriangleIntersection, float distanceToLightPosition, glm::vec3 normal, glm::vec3 lightPosition, glm::vec3 cameraPosition) {

	// ----------------------------- // PROXIMITY

	float lightIntensity =  1.0 / (4.0 * M_PI * (distanceToLightPosition / lightRange) * (distanceToLightPosition / lightRange)); // S / 4*pi*r^2, S = light source strength
	if (lightIntensity > 1.0) {
		lightIntensity = 1.0;
	}
	// if (! proximityOn) {
	// 	lightIntensity = 1;
	// }
	// assert(lightIntensity >= 0.0 && lightIntensity <= 1.0);

	// ----------------------------- // INCIDENCE

	glm::vec3 lightDirection = lightPosition - rayTriangleIntersection.intersectionPoint; // Or lightPosition - newCameraPosition
	lightDirection = glm::normalize(lightDirection);
	float angleOfIncidence = glm::dot(normal, lightDirection);
	if (angleOfIncidence < 0.0) {
		angleOfIncidence = 0.0;
	}
	angleOfIncidence = pow(angleOfIncidence, 0.3); // TODO PLAY HERE

	if (angleOfIncidence > 1.0f) {
		angleOfIncidence = 1.0f;
	}
	// if (! incidenceOn) {
	// 	angleOfIncidence = 1;
	// }
	// assert(angleOfIncidence >= 0.0 && angleOfIncidence <= 1.0);

	// ----------------------------- // SPECULAR

	glm::vec3 viewDirection = cameraPosition - rayTriangleIntersection.intersectionPoint;
	viewDirection = glm::normalize(viewDirection);
	glm::vec3 lightIncidenceDirection = rayTriangleIntersection.intersectionPoint - lightPosition;
	lightIncidenceDirection = glm::normalize(lightIncidenceDirection);
	glm::vec3 lightReflectionDirection = lightIncidenceDirection - (2.0f * normal) * (glm::dot(lightIncidenceDirection, normal));
	glm::normalize(lightReflectionDirection);
	float dotProduct = glm::dot(viewDirection, lightReflectionDirection);
	float specularIntensity = pow(dotProduct, 256);
	if (specularIntensity > 1.0f) {
		specularIntensity = 1.0f;
	}
	// assert(specularIntensity >= 0.0 && specularIntensity <= 1.0);
	if (specularIntensity < 0.5f) {
		specularIntensity = 0.5f;
	}

	// if (! specularOn) {
	// 	specularIntensity = 1;
	// }

	// ----------------------------- //

	return lightIntensity * angleOfIncidence * specularIntensity;
}

std::vector<glm::vec3> lightPoints;

void drawRayTracedScene(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, glm::vec3 cameraPosition, glm::vec3 rayDirection, glm::mat3 cameraOrientation, glm::vec3 lightPosition) {
	window.clearPixels();

	float imageAspectRatio = window.width / (float)window.height; // Assuming width > height 
	// float minDistance = std::numeric_limits<float>::max();
	lightPoints = generateLightPoints(lightPosition);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float brightness = 0.3f;
			if (ambientLight) {
				brightness = 0.3f;
			}
			else {
				brightness = 0.0f;
			}
			rayDirection.x = (2 * ((x + 0.5) / window.width) - 1) * imageAspectRatio; 
			rayDirection.y = (1 - 2 * ((y + 0.5) / window.height));
			rayDirection.z = -2.0;	

			// Make it look bigger 
			// rayDirection.x /= 3.0;
			// rayDirection.y /= 3.0;

			rayDirection = glm::normalize(cameraOrientation * rayDirection); // TODO: Problem when zooming in and out!

			glm::vec3 barycentricCoordinates(0.0, 0.0, 0.0);
			// Fire a ray from the main camera into the scene and get the closest ray/triangle intersection (if it exists)
			RayTriangleIntersection rayTriangleIntersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles, cameraOrientation, false, barycentricCoordinates);

			// We hit something
			if (rayTriangleIntersection.distanceFromCamera != std::numeric_limits<float>::max()) {
				
				float sumBrightness = 0.0f;

				for (int lightIndex = 0; lightIndex < lightPoints.size(); lightIndex++) {

					lightPosition = lightPoints[lightIndex];

					glm::vec3 newCameraPosition = rayTriangleIntersection.intersectionPoint;
					glm::mat3 newCameraOrientation = lookAt(newCameraPosition, lightPosition);
					glm::vec3 newRayDirection = lightPosition - newCameraPosition; // ?
					glm::normalize(newRayDirection); // TODO: Problem when zooming in and out!

					float distanceToLightPosition = glm::distance(lightPosition, newCameraPosition);

					maxDistanceToLightPosition = std::max(maxDistanceToLightPosition, distanceToLightPosition); // 1.05489
					minDistanceToLightPosition = std::min(minDistanceToLightPosition, distanceToLightPosition); // 0.115694

					// Fire a ray from the object to the light source and check if it hits it (i.e. Can the object see the light?)
					RayTriangleIntersection newRayTriangleIntersection = getClosestIntersection(newCameraPosition, newRayDirection, modelTriangles, newCameraOrientation, true, barycentricCoordinates);

					// Light can't be seen
					if (newRayTriangleIntersection.distanceFromCamera < distanceToLightPosition && rayTriangleIntersection.triangleIndex != newRayTriangleIntersection.triangleIndex) { 
						// window.setPixelColour(x, y, 0);
					} 
					else {
						// Light can't reach that surface area
						if (distanceToLightPosition > lightRange) {

						}
						else {	

							if (! gouraud && ! phong) {
								// Calculate normal vector
								rayTriangleIntersection.intersectedTriangle.normal = calculateTriangleNormal(rayTriangleIntersection.intersectedTriangle);
								rayTriangleIntersection.intersectedTriangle.normal = glm::normalize(rayTriangleIntersection.intersectedTriangle.normal);
								brightness = calculateBrightness(rayTriangleIntersection, distanceToLightPosition, rayTriangleIntersection.intersectedTriangle.normal, lightPosition, cameraPosition);
							}
							else if (gouraud) {
								if (rayTriangleIntersection.intersectedTriangle.hasVertexNormals) {
									std::vector<float> brightnesses;
									
									// Calculate brightness for each vertex normal
									for (int i = 0; i < rayTriangleIntersection.intersectedTriangle.vertices.size(); i++) {
										brightnesses.push_back(calculateBrightness(rayTriangleIntersection, distanceToLightPosition, rayTriangleIntersection.intersectedTriangle.verticesNormal[i], lightPosition, cameraPosition));
									}
									// Interpolate brightness with barycentric coordiantes
									brightness = barycentricCoordinates.x * brightnesses[1] + barycentricCoordinates.y * brightnesses[2] + barycentricCoordinates.z * brightnesses[0];
								}
							}
							else if (phong) {
								if (rayTriangleIntersection.intersectedTriangle.hasVertexNormals) {
									glm::vec3 normal = calculateTriangleNormalFromBarycentricCoordinates(rayTriangleIntersection.intersectedTriangle, barycentricCoordinates);
									normal = glm::normalize(normal);
									brightness = calculateBrightness(rayTriangleIntersection, distanceToLightPosition, normal, lightPosition, cameraPosition);
								}
							}

							assert(brightness <= 1.0f);

							if (brightness < 0.3f && ambientLight) {
								brightness = 0.3f;
							}
							sumBrightness += brightness;
						}
					}
				}

				brightness = sumBrightness / (float) lightPoints.size();
				if (brightness < 0.3f && ambientLight) {
					brightness = 0.3f;
				}
				if (brightness > 1.0f) {
					brightness = 1.0f;
				}

				if (defaultLight) {
					brightness = 1.0f;
				}

				Colour RGBColour = rayTriangleIntersection.intersectedTriangle.colour;
				uint32_t colour = (255 << 24) + (int(RGBColour.red * brightness) << 16) + (int(RGBColour.green * brightness) << 8) + int(RGBColour.blue * brightness);

				if (rayTriangleIntersection.intersectedTriangle.hasTexture) {

					float xCoord = rayTriangleIntersection.intersectedTriangle.texturePoints[0].x * barycentricCoordinates.z + 
								rayTriangleIntersection.intersectedTriangle.texturePoints[1].x * barycentricCoordinates.x +
								rayTriangleIntersection.intersectedTriangle.texturePoints[2].x * barycentricCoordinates.y;
					
					float yCoord = rayTriangleIntersection.intersectedTriangle.texturePoints[0].y * barycentricCoordinates.z + 
								rayTriangleIntersection.intersectedTriangle.texturePoints[1].y * barycentricCoordinates.x +
								rayTriangleIntersection.intersectedTriangle.texturePoints[2].y * barycentricCoordinates.y;
					
					uint32_t textureColour;
					if (rayTriangleIntersection.intersectedTriangle.colour.name == "Logo") {
						textureColour = textureMap.pixels[textureMap.width * round(yCoord) + round(xCoord)];
					}
					else if (rayTriangleIntersection.intersectedTriangle.colour.name == "Chess") {
						textureColour = chessTextureMap.pixels[chessTextureMap.width * round(yCoord) + round(xCoord)];
					}
					else if (rayTriangleIntersection.intersectedTriangle.colour.name == "Cobbles") {
						textureColour = cobblesTextureMap.pixels[cobblesTextureMap.width * round(yCoord) + round(xCoord)];
					}
					
					// uint32_t textureColour = textureMap.pixels[textureMap.width * round(yCoord) + round(xCoord)];
					// uint32_t textureColour = textureMap.pixels[textureMap.width * round(yCoord) + round(xCoord)];

					colour = (255 << 24) + (int(((textureColour >> 16) & 0xFF) * brightness) << 16) + (int(((textureColour >> 8) & 0xFF) * brightness) << 8) + int((textureColour & 0xFF) * brightness);
				}

				window.setPixelColour(x, y, colour);
			}
		}
	}
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfSteps) {
	std::vector<glm::vec3> result;
	glm::vec3 ratio = (to - from) /= (numberOfSteps - 1);
	result.push_back(from);
	for (size_t i = 1; i < numberOfSteps; i++) {
		result.push_back(result[i - 1] + ratio);
	}
	return result;
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

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, bool &enableOrbit, glm::vec3 &lightPosition) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_1) { renderOption = Wireframe; }
		else if (event.key.keysym.sym == SDLK_2) { renderOption = Rasterise; }
		else if (event.key.keysym.sym == SDLK_3) { renderOption = Raytrace; }
		else if (event.key.keysym.sym == SDLK_o) { enableOrbit = !enableOrbit; }
		else if (event.key.keysym.sym == SDLK_p) { lightPosition.y += 0.05; }
		else if (event.key.keysym.sym == SDLK_l) { lightPosition.y -= 0.05; }
		else if (event.key.keysym.sym == SDLK_i) { lightPosition.x += 0.05; }
		else if (event.key.keysym.sym == SDLK_j) { lightPosition.x -= 0.05; }
		else if (event.key.keysym.sym == SDLK_u) { lightPosition.z -= 0.05; }
		else if (event.key.keysym.sym == SDLK_h) { lightPosition.z += 0.05; }
		else if (event.key.keysym.sym == SDLK_c) {
			std::cout << "Camera Position: " << "(" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
			std::cout << "Light Position: " << "(" << lightPosition.x << ", " << lightPosition.y << ", " << lightPosition.z << ")" << std::endl;
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
		else if (event.key.keysym.sym == SDLK_b) {
			cameraPosition = rotateCameraAboutYAxis(cameraPosition, 3.14 / 16.0);
			cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
		}	
		else if (event.key.keysym.sym == SDLK_v) {
			cameraPosition = rotateCameraAboutXAxis(cameraPosition, 3.14 / 16.0);
			cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
		}			
	} 
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	// std::string sphereObjFilename = "./models/sphere.obj";
	std::string objFilename = "./models/cornell-box.obj";
	// std::string objFilename = "./models/textured-cornell-box.obj";
	// std::string objFilename = "./models/logo.obj";
	// std::string objFilename = "./models/modified2.obj";
	// std::string objFilename = sphereObjFilename;

	std::string mtlFilename = "./models/cornell-box.mtl";
	// std::string mtlFilename = "./models/textured-cornell-box.mtl";

	std::vector<Colour> colourPalette;
	std::vector<ModelTriangle> modelTriangles;
	std::vector<std::vector<float> > depthBuffer(WIDTH, std::vector<float> (HEIGHT, 0));

	readMTLFile(mtlFilename, colourPalette);
	parseOBJFile(objFilename, colourPalette, modelTriangles);

	// for (int i = 0; i < modelTriangles.size(); i++) {
	// 	std::cout << i << " " << modelTriangles[i].hasTexture << std::endl;
	// 		// std::cout << "[0]" << "(" << modelTriangles[i].verticesNormal[0].x  << ", " << modelTriangles[i].verticesNormal[0].y << ", " << modelTriangles[i].verticesNormal[0].z << ")" << std::endl;
	// 		// std::cout << "[1]" << "(" << modelTriangles[i].verticesNormal[1].x  << ", " << modelTriangles[i].verticesNormal[1].y << ", " << modelTriangles[i].verticesNormal[1].z << ")" << std::endl;
	// 		// std::cout << "[2]" << "(" << modelTriangles[i].verticesNormal[2].x  << ", " << modelTriangles[i].verticesNormal[2].y << ", " << modelTriangles[i].verticesNormal[2].z << ")" << std::endl;
	// 		// std::cout << "-------------------------------" << std::endl;
	// }

	glm::vec3 cameraPosition(0.0, 0.0, 2.0); // Try (0.0, 0.0, 4.0) for rasterising
	glm::vec3 rayDirection(0.0, 0.0, -1.0);
	glm::vec3 lightPosition(0.0, 0.35, 0.15); // For cornell box
	// glm::vec3 lightPosition(-0.45, 0.65, 0.85); // For spere

	lightPoints = generateLightPoints(lightPosition);

	glm::vec3 right(1.0, 0.0, 0.0); 
	glm::vec3 up(0.0, 1.0, 0.0);
	glm::vec3 forward(0.0, 0.0, 1.0);

	glm::mat3 cameraOrientation(right, up, forward);

	// lightPosition.y = 0.30;
	

	// glm::vec3 column0 = glm::vec3(1.0, 0.0, 0.0);
	// glm::vec3 column1 = glm::vec3(0.0, 0.7, 0.7);
	// glm::vec3 column2 = glm::vec3(0.0, -0.7, 0.7);
	// glm::mat3 n = glm::mat3(column0, column1, column2); 

	// glm::mat3 rotationDefault = glm::mat3(
	// 	1.0, 0.0, 0.0,
	// 	0.0, 0.0, 1.0,
	// 	0.0, 1.0, 0.0
	// );

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

	int frame = 0;

	// std::vector<glm::vec3> cameraPositions = interpolateThreeElementValues(cameraPosition, glm::vec3(0.0, 0.0, 0.0), 100);
	// for (int i = 0; i < 100; i++) {
	// 	cameraPosition = cameraPositions[i];
	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	window.renderFrame();
	// 	std::string filename = "animation1/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
	// }
	// return 0;

	// std::vector<glm::vec3> cameraPositions = interpolateThreeElementValues(cameraPosition, glm::vec3(0.0, 0.0, 0.0), 100);
	// for (int i = 0; i < 190; i++) {
	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
		
	// 	// float angle = 3.14 / 16.0;
	// 	if (i < 50) {
	// 		drawWireframedScene(window, modelTriangles, cameraPosition, cameraOrientation, depthBuffer);
	// 		float angle = 3.14 / 16.0;
	// 		cameraPosition = rotateCameraAboutYAxis(cameraPosition, angle);
	// 		cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	}
	// 	else if (i >= 50 && i < 100) {
	// 		drawRasterizedScene(window, modelTriangles, cameraPosition, cameraOrientation, depthBuffer);
	// 		float angle = 3.14 / 16.0;
	// 		cameraPosition = rotateCameraAboutYAxis(cameraPosition, angle);
	// 		cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	}
	// 	else {
	// 		drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 		float angle = 3.14 / 16.0;
	// 		cameraPosition = rotateCameraAboutYAxis(cameraPosition, angle);
	// 		cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	}

	// 	window.renderFrame();
	// 	std::string filename = "animation2/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
	// }
	// return 0;

	// std::vector<glm::vec3> cameraPositions = interpolateThreeElementValues(cameraPosition, glm::vec3(0.0, 0.0, 0.2), 100);
	// for (int i = 0; i < 100; i++) {
	// 	cameraPosition = cameraPositions[i];
	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	window.renderFrame();
	// 	std::string filename = "animation3/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
	// }
	// return 0;
	
	// std::vector<glm::vec3> cameraPositions = interpolateThreeElementValues(glm::vec3(0.0, 0.0, 0.2), glm::vec3(2.0, 0.0, 0.0), 50);
	// for (int i = 100; i < 150; i++) {
	// 	cameraPosition = cameraPositions[i];
	// 	cameraOrientation = lookAt(cameraPosition, glm::vec3(2.0, 0.0, 0.0));
	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	window.renderFrame();
	// 	std::string filename = "animation3/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
		
	// }
	// return 0;

	// std::vector<glm::vec3> cameraPositions = interpolateThreeElementValues(glm::vec3(0.0, 0.0, -0.4), glm::vec3(0.0, 0.0, 2.0), 150);
	// for (int i = 0; i < 150; i++) {
	// 	cameraPosition = cameraPositions[i];
	// 	// cameraOrientation = lookAt(cameraPosition, glm::vec3(2.0, 0.0, 0.0));
	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	window.renderFrame();
	// 	std::string filename = "animation4/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
		
	// }
	// return 0;

	// for (int i = 0; i < 20; i++) {
	// 	// cameraPosition = cameraPositions[i];
	// 	// cameraOrientation = lookAt(cameraPosition, glm::vec3(2.0, 0.0, 0.0));
	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);		
	// 	window.renderFrame();
	// 	std::string filename = "animation7/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
		
	// }
	// return 0;

	// for (int i = 0; i < 22; i++) {
	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
		
	// 	// float angle = 3.14 / 16.0;

	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	float angle = 3.14 / 128.0;
	// 	cameraPosition = rotateCameraAboutYAxis(cameraPosition, angle);
	// 	cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	// window.renderFrame();
	// 	// std::string filename = "animation8/" + std::to_string(i) + ".ppm";
	// 	// window.savePPM(filename);
	// }

	// std::cout << cameraPosition.x << cameraPosition.y << cameraPosition. z << std::endl;

	// for (int i = 21; i < 39; i++) {
	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
		
	// 	// float angle = 3.14 / 16.0;

	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	float angle = 3.14 / 128.0;
	// 	cameraPosition = rotateCameraAboutYAxis(cameraPosition, -angle);
	// 	cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	// window.renderFrame();
	// 	// // std::string filename = "animation8/" + std::to_string(i) + ".ppm";
	// 	// window.savePPM(filename);
	// }
	// 	std::cout << cameraPosition.x << cameraPosition.y << cameraPosition. z << std::endl;


	// for (int i = 39; i < 21 * 3; i++) {
	// 	// drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
		
	// 	// float angle = 3.14 / 16.0;

	// 	drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition);
	// 	float angle = 3.14 / 128.0;
	// 	cameraPosition = rotateCameraAboutYAxis(cameraPosition, -angle);
	// 	cameraOrientation = lookAt(cameraPosition, glm::vec3(0.0, 0.0, 0.0));
	// 	window.renderFrame();
	// 	std::string filename = "animation8/" + std::to_string(i) + ".ppm";
	// 	window.savePPM(filename);
	// }
	// return 0;

	// cameraPosition = glm::vec3(0.0, 0.0, 0.0);
	// renderOption = Raytrace;
	

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cameraPosition, cameraOrientation, enableOrbit, lightPosition);
		switch (renderOption) {
			case Wireframe: drawWireframedScene(window, modelTriangles, cameraPosition, cameraOrientation, depthBuffer); break;
			case Rasterise: drawRasterizedScene(window, modelTriangles, cameraPosition, cameraOrientation, depthBuffer); break;
			case Raytrace: drawRayTracedScene(window, modelTriangles, cameraPosition, rayDirection, cameraOrientation, lightPosition); break;
		}

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
		// frame++;
	}
}
