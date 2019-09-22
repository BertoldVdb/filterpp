#ifndef KEEPER_H_
#define KEEPER_H_

#include <cstring>
#include <vector>

namespace Util {

template <typename T> class Keeper {
public:
	Keeper(size_t bufferSize):
			buffer_(bufferSize){
	}

	inline void insert(T* samples, size_t numSamples){
		size_t index = 0;
		auto len = numSamples;

		while(len){
			auto write = len;
			if(write > buffer_.size()-bufferIndex_){
				write = buffer_.size()-bufferIndex_;
			}

			std::memcpy(&buffer_[bufferIndex_], samples+index, write*sizeof(T));
			len -= write;
			index += write;
			bufferIndex_ += write;

			if(bufferIndex_>=buffer_.size()){
				bufferIndex_ -= buffer_.size();
			}
		}

		bufferLevel_ += numSamples;
		if(bufferLevel_ > buffer_.size()){
			bufferLevel_ = buffer_.size();
		}
	}

	inline void extract(T* samples){
		std::memcpy(samples, &buffer_[bufferIndex_], (buffer_.size() - bufferIndex_)*sizeof(T));
		std::memcpy(samples + buffer_.size() - bufferIndex_, &buffer_[0], bufferIndex_*sizeof(T));
	}

	inline size_t size(){
		return buffer_.size();
	}

	inline size_t remaining(){
		return buffer_.size() - bufferLevel_;
	}

	inline size_t level(){
		return bufferLevel_;
	}

	inline void clear(){
		bufferLevel_ = 0;
	}

	inline void drop(size_t samples){
		if(samples >= bufferLevel_){
			bufferLevel_ = 0;
			return;
		}

		bufferLevel_ -= samples;
	}

private:
	size_t bufferIndex_ = 0;
	size_t bufferLevel_ = 0;

	std::vector<T> buffer_;
};

}

#endif /* KEEPER_H_ */
