import Heart from './svg/heart.svg';

interface IMouseIcon {
  size: number;
}

function HeartIcon({ size }: IMouseIcon) {
  return (
    <div style={{ width: size }}>
      <img src={Heart} alt="heart icon" />
    </div>
  );
}

export default HeartIcon;
